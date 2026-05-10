#include "ccalc.h"
#include "imgui.h"
#include "imgui_impl_win32.h"
#include "imgui_impl_dx11.h"
#include <d3d11.h>
#include <tchar.h>
#include <cmath>
#include <vector>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

static ID3D11Device*            g_pd3dDevice = nullptr;
static ID3D11DeviceContext*     g_pd3dDeviceContext = nullptr;
static IDXGISwapChain*          g_pSwapChain = nullptr;
static bool                     g_SwapChainOccluded = false;
static UINT                     g_ResizeWidth = 0, g_ResizeHeight = 0;
static ID3D11RenderTargetView*  g_mainRenderTargetView = nullptr;

bool CreateDeviceD3D(HWND hWnd);
void CleanupDeviceD3D();
void CreateRenderTarget();
void CleanupRenderTarget();
LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

enum GraphType { GT_EXPLICIT, GT_IMPLICIT, GT_PARAMETRIC, GT_POLAR };


struct EvalVars {
    double x = 0, y = 0, t = 0, theta = 0;
};

struct FuncEntry {
    GraphType gtype;
    char expr[256];
    char expr2[256];
    double t_min, t_max;
    bool visible;
    ImU32 color;
    bool valid;
    std::string error_msg;
    ASTPtr cached_ast;
    ASTPtr cached_ast2;
    std::string cached_expr_str;
    std::string cached_expr2_str;
};

static const ImU32 func_colors[] = {
    IM_COL32(255, 100, 100, 255),
    IM_COL32(100, 255, 100, 255),
    IM_COL32(100, 150, 255, 255),
    IM_COL32(255, 255, 100, 255),
    IM_COL32(255, 100, 255, 255),
    IM_COL32(100, 255, 255, 255),
    IM_COL32(255, 180, 100, 255),
    IM_COL32(180, 100, 255, 255),
};

static bool ensure_ast(FuncEntry& f) {
    bool need_reparse = false;
    std::string cur(f.expr);
    if (cur != f.cached_expr_str) need_reparse = true;
    if (f.gtype == GT_PARAMETRIC) {
        std::string cur2(f.expr2);
        if (cur2 != f.cached_expr2_str) need_reparse = true;
    }
    if (!need_reparse && f.cached_ast) return f.valid;

    f.cached_expr_str = cur;
    f.cached_expr2_str = f.expr2;
    f.cached_ast.reset();
    f.cached_ast2.reset();
    f.valid = true;
    f.error_msg = "";

    if (cur.empty()) return true;

    try {
        Lexer lexer(cur);
        auto tokens = lexer.tokenize();
        Parser parser(tokens);
        f.cached_ast = parser.parse();
    } catch (const std::exception& e) {
        f.valid = false;
        f.error_msg = std::string("expr1: ") + e.what();
        return false;
    }

    if (f.gtype == GT_PARAMETRIC && f.expr2[0] != '\0') {
        try {
            Lexer lexer2(f.expr2);
            auto tokens2 = lexer2.tokenize();
            Parser parser2(tokens2);
            f.cached_ast2 = parser2.parse();
        } catch (const std::exception& e) {
            f.valid = false;
            f.error_msg = std::string("expr2: ") + e.what();
            return false;
        }
    }

    return true;
}

static double fast_eval(ASTPtr node, const EvalVars& v, bool deg_mode) {
    if (!node) return NAN;
    switch (node->type) {
    case ASTNode::NUMBER:
        return std::stod(node->number.to_string());
    case ASTNode::CONSTANT:
        if (node->name == "pi") return M_PI;
        if (node->name == "e") return M_E;
        return NAN;
    case ASTNode::VARIABLE:
        if (node->name == "x") return v.x;
        if (node->name == "y") return v.y;
        if (node->name == "t") return v.t;
        if (node->name == "theta") return v.theta;
        return 0;
    case ASTNode::BINOP: {
        double l = fast_eval(node->left, v, deg_mode);
        double r = fast_eval(node->right, v, deg_mode);
        switch (node->op) {
        case '+': return l + r;
        case '-': return l - r;
        case '*': return l * r;
        case '/': return (r == 0) ? NAN : l / r;
        case '^': return pow(l, r);
        case '%': return (r == 0) ? NAN : fmod(l, r);
        }
        return NAN;
    }
    case ASTNode::UNARYOP:
        if (node->op == '-') return -fast_eval(node->left, v, deg_mode);
        return fast_eval(node->left, v, deg_mode);
    case ASTNode::FACTORIAL: {
        double val = fast_eval(node->left, v, deg_mode);
        if (val < 0 || val != floor(val) || val > 170) return NAN;
        double r = 1;
        for (int i = 2; i <= (int)val; i++) r *= i;
        return r;
    }
    case ASTNode::FUNCTION: {
        auto& args = node->args;
        auto fe = [&](int i) { return fast_eval(args[i], v, deg_mode); };
        const std::string& name = node->name;
        auto to_rad = [&](double val) { return deg_mode ? val * M_PI / 180.0 : val; };
        auto from_rad = [&](double val) { return deg_mode ? val * 180.0 / M_PI : val; };

        if (name == "sin" && args.size() == 1) return sin(to_rad(fe(0)));
        if (name == "cos" && args.size() == 1) return cos(to_rad(fe(0)));
        if (name == "tan" && args.size() == 1) { double c = cos(to_rad(fe(0))); return fabs(c) < 1e-15 ? NAN : tan(to_rad(fe(0))); }
        if (name == "asin" && args.size() == 1) return from_rad(asin(fe(0)));
        if (name == "acos" && args.size() == 1) return from_rad(acos(fe(0)));
        if (name == "atan" && args.size() == 1) return from_rad(atan(fe(0)));
        if (name == "sinh" && args.size() == 1) return sinh(fe(0));
        if (name == "cosh" && args.size() == 1) return cosh(fe(0));
        if (name == "tanh" && args.size() == 1) return tanh(fe(0));
        if (name == "sqrt" && args.size() == 1) return fe(0) < 0 ? NAN : sqrt(fe(0));
        if (name == "cbrt" && args.size() == 1) return cbrt(fe(0));
        if (name == "nrt" && args.size() == 2) { double n = fe(0); return n == 0 ? NAN : pow(fe(1), 1.0 / n); }
        if (name == "abs" && args.size() == 1) return fabs(fe(0));
        if (name == "ln" && args.size() == 1) return fe(0) <= 0 ? NAN : log(fe(0));
        if (name == "lg" && args.size() == 1) return fe(0) <= 0 ? NAN : log10(fe(0));
        if (name == "log" && args.size() == 1) return fe(0) <= 0 ? NAN : log(fe(0));
        if (name == "log" && args.size() == 2) { double b = fe(0); double val = fe(1); return (b <= 0 || val <= 0) ? NAN : log(val) / log(b); }
        if (name == "exp" && args.size() == 1) return exp(fe(0));
        if (name == "fact" && args.size() == 1) { double val = fe(0); if (val < 0 || val != floor(val) || val > 170) return NAN; double r = 1; for (int i = 2; i <= (int)val; i++) r *= i; return r; }
        if (name == "floor" && args.size() == 1) return floor(fe(0));
        if (name == "ceil" && args.size() == 1) return ceil(fe(0));
        if (name == "round" && args.size() == 1) return round(fe(0));
        if (name == "sign" && args.size() == 1) { double val = fe(0); return val > 0 ? 1 : (val < 0 ? -1 : 0); }
        if (name == "max" && args.size() >= 2) { double r = fe(0); for (size_t i = 1; i < args.size(); i++) r = fmax(r, fe(i)); return r; }
        if (name == "min" && args.size() >= 2) { double r = fe(0); for (size_t i = 1; i < args.size(); i++) r = fmin(r, fe(i)); return r; }
        if (name == "pow" && args.size() == 2) return pow(fe(0), fe(1));
        if (name == "root" && args.size() == 2) { double n = fe(0); return n == 0 ? NAN : pow(fe(1), 1.0 / n); }
        if (name == "log2" && args.size() == 1) return fe(0) <= 0 ? NAN : log2(fe(0));
        if (name == "hypot" && args.size() == 2) return hypot(fe(0), fe(1));
        if (name == "atan2" && args.size() == 2) return from_rad(atan2(fe(0), fe(1)));
        if (name == "deg" && args.size() == 1) return fe(0) * M_PI / 180.0;
        if (name == "rad" && args.size() == 1) return fe(0) * 180.0 / M_PI;
        return NAN;
    }
    case ASTNode::VEC_LITERAL:
        return NAN;
    }
    return NAN;
}

static double nice_step(double range) {
    if (range <= 0) return 1;
    double rough = range / 10.0;
    double mag = pow(10.0, floor(log10(rough)));
    double norm = rough / mag;
    double nice;
    if (norm < 1.5) nice = 1;
    else if (norm < 3.0) nice = 2;
    else if (norm < 7.0) nice = 5;
    else nice = 10;
    return nice * mag;
}

static std::string format_tick(double val) {
    if (fabs(val) < 1e-15) return "0";
    if (fabs(val) >= 1e6 || (fabs(val) < 0.001 && fabs(val) > 0)) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.1e", val);
        return buf;
    }
    char buf[32];
    snprintf(buf, sizeof(buf), "%g", val);
    return buf;
}

static void draw_explicit(ImDrawList* dl, FuncEntry& f, bool deg_mode,
    double x_min, double x_max, double y_min, double y_max,
    const std::function<ImVec2(double, double)>& w2s, float canvas_w) {
    if (!ensure_ast(f) || !f.cached_ast) return;
    int n_samples = (int)canvas_w;
    if (n_samples < 100) n_samples = 100;
    if (n_samples > 2000) n_samples = 2000;
    double dx = (x_max - x_min) / n_samples;
    std::vector<ImVec2> points;
    points.reserve(n_samples + 1);
    double prev_y = NAN;
    for (int i = 0; i <= n_samples; i++) {
        double xv = x_min + i * dx;
        EvalVars v; v.x = xv;
        double yv = fast_eval(f.cached_ast, v, deg_mode);
        if (!std::isfinite(yv) || fabs(yv) > 1e10) {
            if (points.size() >= 2) dl->AddPolyline(points.data(), (int)points.size(), f.color, false, 2.0f);
            points.clear(); prev_y = NAN; continue;
        }
        if (std::isfinite(prev_y) && fabs(yv - prev_y) > (y_max - y_min) * 10) {
            if (points.size() >= 2) dl->AddPolyline(points.data(), (int)points.size(), f.color, false, 2.0f);
            points.clear();
        }
        points.push_back(w2s(xv, yv));
        prev_y = yv;
    }
    if (points.size() >= 2) dl->AddPolyline(points.data(), (int)points.size(), f.color, false, 2.0f);
}

static void draw_implicit(ImDrawList* dl, FuncEntry& f, bool deg_mode,
    double x_min, double x_max, double y_min, double y_max,
    const std::function<ImVec2(double, double)>& w2s, float canvas_w, float canvas_h) {
    if (!ensure_ast(f) || !f.cached_ast) return;
    int nx = std::min(300, std::max(100, (int)canvas_w / 3));
    int ny = std::min(300, std::max(100, (int)canvas_h / 3));
    double dx = (x_max - x_min) / nx;
    double dy = (y_max - y_min) / ny;
    std::vector<std::vector<double>> grid(ny + 1, std::vector<double>(nx + 1));
    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i <= nx; i++) {
            EvalVars v;
            v.x = x_min + i * dx;
            v.y = y_min + j * dy;
            grid[j][i] = fast_eval(f.cached_ast, v, deg_mode);
        }
    }
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double f00 = grid[j][i], f10 = grid[j][i + 1];
            double f01 = grid[j + 1][i], f11 = grid[j + 1][i + 1];
            int code = 0;
            if (f00 > 0) code |= 1;
            if (f10 > 0) code |= 2;
            if (f11 > 0) code |= 4;
            if (f01 > 0) code |= 8;
            if (code == 0 || code == 15) continue;
            double x0 = x_min + i * dx, x1 = x0 + dx;
            double y0 = y_min + j * dy, y1 = y0 + dy;
            auto lerp = [](double a, double fa, double b, double fb) -> double {
                if (fabs(fa - fb) < 1e-15) return (a + b) / 2;
                return a + (b - a) * fa / (fa - fb);
            };
            struct EP { double x, y; };
            EP edges[4] = {};
            bool he[4] = {};
            if ((f00 > 0) != (f10 > 0)) { edges[0] = {lerp(x0, f00, x1, f10), y0}; he[0] = true; }
            if ((f10 > 0) != (f11 > 0)) { edges[1] = {x1, lerp(y0, f10, y1, f11)}; he[1] = true; }
            if ((f01 > 0) != (f11 > 0)) { edges[2] = {lerp(x0, f01, x1, f11), y1}; he[2] = true; }
            if ((f00 > 0) != (f01 > 0)) { edges[3] = {x0, lerp(y0, f00, y1, f01)}; he[3] = true; }
            int cnt = he[0] + he[1] + he[2] + he[3];
            if (cnt == 2) {
                EP a = {}, b = {};
                bool first = true;
                for (int k = 0; k < 4; k++) {
                    if (he[k]) {
                        if (first) { a = edges[k]; first = false; }
                        else { b = edges[k]; }
                    }
                }
                dl->AddLine(w2s(a.x, a.y), w2s(b.x, b.y), f.color, 2.0f);
            } else if (cnt == 4) {
                double avg = (f00 + f10 + f01 + f11) / 4.0;
                if ((avg > 0) == (f00 > 0)) {
                    dl->AddLine(w2s(edges[0].x, edges[0].y), w2s(edges[3].x, edges[3].y), f.color, 2.0f);
                    dl->AddLine(w2s(edges[1].x, edges[1].y), w2s(edges[2].x, edges[2].y), f.color, 2.0f);
                } else {
                    dl->AddLine(w2s(edges[0].x, edges[0].y), w2s(edges[1].x, edges[1].y), f.color, 2.0f);
                    dl->AddLine(w2s(edges[2].x, edges[2].y), w2s(edges[3].x, edges[3].y), f.color, 2.0f);
                }
            }
        }
    }
}

static void draw_parametric(ImDrawList* dl, FuncEntry& f, bool deg_mode,
    double, double, double, double,
    const std::function<ImVec2(double, double)>& w2s) {
    if (!ensure_ast(f) || !f.cached_ast || !f.cached_ast2) return;
    int n_samples = 2000;
    double dt = (f.t_max - f.t_min) / n_samples;
    std::vector<ImVec2> points;
    points.reserve(n_samples + 1);
    for (int i = 0; i <= n_samples; i++) {
        double tv = f.t_min + i * dt;
        EvalVars v; v.t = tv;
        double xv = fast_eval(f.cached_ast, v, deg_mode);
        double yv = fast_eval(f.cached_ast2, v, deg_mode);
        if (!std::isfinite(xv) || !std::isfinite(yv) || fabs(xv) > 1e10 || fabs(yv) > 1e10) {
            if (points.size() >= 2) dl->AddPolyline(points.data(), (int)points.size(), f.color, false, 2.0f);
            points.clear(); continue;
        }
        ImVec2 sp = w2s(xv, yv);
        if (!points.empty()) {
            float ddx = sp.x - points.back().x, ddy = sp.y - points.back().y;
            if (ddx * ddx + ddy * ddy > (800 * 800)) {
                if (points.size() >= 2) dl->AddPolyline(points.data(), (int)points.size(), f.color, false, 2.0f);
                points.clear();
            }
        }
        points.push_back(sp);
    }
    if (points.size() >= 2) dl->AddPolyline(points.data(), (int)points.size(), f.color, false, 2.0f);
}

static void draw_polar(ImDrawList* dl, FuncEntry& f, bool deg_mode,
    double, double, double, double,
    const std::function<ImVec2(double, double)>& w2s) {
    if (!ensure_ast(f) || !f.cached_ast) return;
    int n_samples = 2000;
    double dtheta = (f.t_max - f.t_min) / n_samples;
    std::vector<ImVec2> points;
    points.reserve(n_samples + 1);
    for (int i = 0; i <= n_samples; i++) {
        double th = f.t_min + i * dtheta;
        EvalVars v; v.theta = th;
        double r = fast_eval(f.cached_ast, v, deg_mode);
        if (!std::isfinite(r) || fabs(r) > 1e10) {
            if (points.size() >= 2) dl->AddPolyline(points.data(), (int)points.size(), f.color, false, 2.0f);
            points.clear(); continue;
        }
        double xv = r * cos(th);
        double yv = r * sin(th);
        ImVec2 sp = w2s(xv, yv);
        if (!points.empty()) {
            float ddx = sp.x - points.back().x, ddy = sp.y - points.back().y;
            if (ddx * ddx + ddy * ddy > (800 * 800)) {
                if (points.size() >= 2) dl->AddPolyline(points.data(), (int)points.size(), f.color, false, 2.0f);
                points.clear();
            }
        }
        points.push_back(sp);
    }
    if (points.size() >= 2) dl->AddPolyline(points.data(), (int)points.size(), f.color, false, 2.0f);
}

int main(int argc, char** argv) {
    std::vector<FuncEntry> init_funcs;
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "-i" && i + 1 < argc) {
            FuncEntry fe = {}; fe.gtype = GT_IMPLICIT;
            strncpy(fe.expr, argv[++i], 255); fe.t_min = 0; fe.t_max = 2 * M_PI;
            fe.visible = true; fe.color = func_colors[init_funcs.size() % 8];
            init_funcs.push_back(fe);
        } else if (arg == "-p" && i + 2 < argc) {
            FuncEntry fe = {}; fe.gtype = GT_PARAMETRIC;
            strncpy(fe.expr, argv[++i], 255);
            strncpy(fe.expr2, argv[++i], 255);
            fe.t_min = 0; fe.t_max = 2 * M_PI;
            fe.visible = true; fe.color = func_colors[init_funcs.size() % 8];
            init_funcs.push_back(fe);
        } else if (arg == "-l" && i + 1 < argc) {
            FuncEntry fe = {}; fe.gtype = GT_POLAR;
            strncpy(fe.expr, argv[++i], 255); fe.t_min = 0; fe.t_max = 2 * M_PI;
            fe.visible = true; fe.color = func_colors[init_funcs.size() % 8];
            init_funcs.push_back(fe);
        } else if (arg[0] != '-') {
            FuncEntry fe = {}; fe.gtype = GT_EXPLICIT;
            strncpy(fe.expr, argv[i], 255); fe.t_min = 0; fe.t_max = 2 * M_PI;
            fe.visible = true; fe.color = func_colors[init_funcs.size() % 8];
            init_funcs.push_back(fe);
        }
    }

    ImGui_ImplWin32_EnableDpiAwareness();
    float main_scale = ImGui_ImplWin32_GetDpiScaleForMonitor(
        ::MonitorFromPoint(POINT{0, 0}, MONITOR_DEFAULTTOPRIMARY));

    WNDCLASSEXW wc = {sizeof(wc), CS_CLASSDC, WndProc, 0L, 0L,
        GetModuleHandle(nullptr), nullptr, nullptr, nullptr, nullptr,
        L"CCalc Graph", nullptr};
    ::RegisterClassExW(&wc);

    std::wstring wtitle = L"CCalc Graph";
    if (!init_funcs.empty()) {
        std::string title_str;
        for (auto& f : init_funcs) {
            if (!title_str.empty()) title_str += ", ";
            if (f.gtype == GT_EXPLICIT) title_str += "y=" + std::string(f.expr);
            else if (f.gtype == GT_IMPLICIT) title_str += std::string(f.expr) + "=0";
            else if (f.gtype == GT_PARAMETRIC) title_str += "(" + std::string(f.expr) + ", " + std::string(f.expr2) + ")";
            else if (f.gtype == GT_POLAR) title_str += "r=" + std::string(f.expr);
        }
        int len = MultiByteToWideChar(CP_UTF8, 0, title_str.c_str(), -1, nullptr, 0);
        std::wstring wexpr(len, 0);
        MultiByteToWideChar(CP_UTF8, 0, title_str.c_str(), -1, &wexpr[0], len);
        wtitle = wexpr;
    }

    HWND hwnd = ::CreateWindowW(wc.lpszClassName, wtitle.c_str(),
        WS_OVERLAPPEDWINDOW, 100, 100,
        (int)(800 * main_scale), (int)(600 * main_scale),
        nullptr, nullptr, wc.hInstance, nullptr);

    if (!CreateDeviceD3D(hwnd)) {
        CleanupDeviceD3D();
        ::UnregisterClassW(wc.lpszClassName, wc.hInstance);
        return 1;
    }

    ::ShowWindow(hwnd, SW_SHOWDEFAULT);
    ::UpdateWindow(hwnd);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.ScaleAllSizes(main_scale);
    style.FontScaleDpi = main_scale;

    ImGui_ImplWin32_Init(hwnd);
    ImGui_ImplDX11_Init(g_pd3dDevice, g_pd3dDeviceContext);

    ImFontConfig font_cfg;
    font_cfg.OversampleH = 2;
    font_cfg.OversampleV = 1;
    io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf",
        16.0f * main_scale, &font_cfg,
        io.Fonts->GetGlyphRangesDefault());
    ImFont* font_cjk = io.Fonts->AddFontFromFileTTF(
        "c:\\Windows\\Fonts\\msyh.ttc",
        16.0f * main_scale, &font_cfg,
        io.Fonts->GetGlyphRangesChineseFull());
    (void)font_cjk;

    bool deg_mode = false;

    std::vector<FuncEntry> funcs;
    if (!init_funcs.empty()) {
        funcs = init_funcs;
    } else {
        funcs.push_back({GT_EXPLICIT, "sin(x)", "", 0, 2 * M_PI, true, func_colors[0], true, "", nullptr, nullptr, "", ""});
    }
    for (int i = (int)funcs.size(); i < 8; i++) {
        funcs.push_back({GT_EXPLICIT, "", "", 0, 2 * M_PI, true, func_colors[i % 8], true, "", nullptr, nullptr, "", ""});
    }

    double x_min = -10.0, x_max = 10.0;
    double y_min = -10.0, y_max = 10.0;
    bool auto_y = true;
    bool show_grid = true;
    bool show_coords = true;
    bool show_settings = false;
    bool dragging = false;
    ImVec2 drag_start;
    double drag_xmin, drag_xmax, drag_ymin, drag_ymax;
    ImVec4 clear_color = ImVec4(0.12f, 0.12f, 0.14f, 1.00f);

    bool done = false;
    while (!done) {
        MSG msg;
        while (::PeekMessage(&msg, nullptr, 0U, 0U, PM_REMOVE)) {
            ::TranslateMessage(&msg);
            ::DispatchMessage(&msg);
            if (msg.message == WM_QUIT) done = true;
        }
        if (done) break;

        if (g_SwapChainOccluded &&
            g_pSwapChain->Present(0, DXGI_PRESENT_TEST) == DXGI_STATUS_OCCLUDED) {
            ::Sleep(10); continue;
        }
        g_SwapChainOccluded = false;

        if (g_ResizeWidth != 0 && g_ResizeHeight != 0) {
            CleanupRenderTarget();
            g_pSwapChain->ResizeBuffers(0, g_ResizeWidth, g_ResizeHeight, DXGI_FORMAT_UNKNOWN, 0);
            g_ResizeWidth = g_ResizeHeight = 0;
            CreateRenderTarget();
        }

        ImGui_ImplDX11_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();

        {
            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::SetNextWindowSize(io.DisplaySize);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
            ImGui::Begin("CCalc Graph", nullptr,
                ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_NoBringToFrontOnFocus);

            if (show_settings) {
                float panel_width = 310.0f * main_scale;
                ImGui::BeginChild("SettingsPanel", ImVec2(panel_width, 0), true);

                ImGui::Text("Functions:");
                for (int i = 0; i < (int)funcs.size(); i++) {
                    ImGui::PushID(i);
                    ImGui::Checkbox("##vis", &funcs[i].visible);
                    ImGui::SameLine();
                    ImGui::ColorButton("##col",
                        ImVec4(((funcs[i].color >> 0) & 0xFF) / 255.0f,
                               ((funcs[i].color >> 8) & 0xFF) / 255.0f,
                               ((funcs[i].color >> 16) & 0xFF) / 255.0f, 1.0f),
                        ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoPicker,
                        ImVec2(16, 16));
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(70);
                    ImGui::Combo("##type", (int*)&funcs[i].gtype, "y=f(x)\0f(x,y)=0\0param\0polar\0");
                    ImGui::PushItemWidth(-1);
                    if (funcs[i].gtype == GT_PARAMETRIC) {
                        ImGui::Text("  x(t):");
                        ImGui::SameLine();
                        ImGui::SetNextItemWidth(-1);
                        ImGui::InputText("##expr", funcs[i].expr, sizeof(funcs[i].expr));
                        ImGui::Text("  y(t):");
                        ImGui::SameLine();
                        ImGui::SetNextItemWidth(-1);
                        ImGui::InputText("##expr2", funcs[i].expr2, sizeof(funcs[i].expr2));
                    } else if (funcs[i].gtype == GT_IMPLICIT) {
                        ImGui::Text("  f(x,y)=0:");
                        ImGui::SameLine();
                        ImGui::SetNextItemWidth(-1);
                        ImGui::InputText("##expr", funcs[i].expr, sizeof(funcs[i].expr));
                    } else if (funcs[i].gtype == GT_POLAR) {
                        ImGui::Text("  r(theta):");
                        ImGui::SameLine();
                        ImGui::SetNextItemWidth(-1);
                        ImGui::InputText("##expr", funcs[i].expr, sizeof(funcs[i].expr));
                    } else {
                        ImGui::InputText("##expr", funcs[i].expr, sizeof(funcs[i].expr));
                    }
                    ImGui::PopItemWidth();
                    if (funcs[i].gtype == GT_PARAMETRIC || funcs[i].gtype == GT_POLAR) {
                        ImGui::SetNextItemWidth(80);
                        ImGui::InputDouble("t min", &funcs[i].t_min, 0, 0, "%.2f");
                        ImGui::SameLine();
                        ImGui::SetNextItemWidth(80);
                        ImGui::InputDouble("t max", &funcs[i].t_max, 0, 0, "%.2f");
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("X##del")) {
                        funcs.erase(funcs.begin() + i); i--;
                        ImGui::PopID(); continue;
                    }
                    ensure_ast(funcs[i]);
                    if (!funcs[i].valid && funcs[i].expr[0] != '\0') {
                        ImGui::TextColored(ImVec4(1, 0.5f, 0.5f, 1), "  %s", funcs[i].error_msg.c_str());
                    }
                    ImGui::PopID();
                }
                if (ImGui::Button("+ Add")) {
                    int ci = (int)funcs.size() % 8;
                    funcs.push_back({GT_EXPLICIT, "", "", 0, 2 * M_PI, true, func_colors[ci], true, "", nullptr, nullptr, "", ""});
                }

                ImGui::Separator();
                ImGui::Checkbox("Auto Y range", &auto_y);
                float fx_min = (float)x_min, fx_max = (float)x_max;
                float fy_min = (float)y_min, fy_max = (float)y_max;
                ImGui::SetNextItemWidth(100);
                if (ImGui::InputFloat("X min", &fx_min, 0, 0, "%.2f")) x_min = fx_min;
                ImGui::SameLine();
                ImGui::SetNextItemWidth(100);
                if (ImGui::InputFloat("X max", &fx_max, 0, 0, "%.2f")) x_max = fx_max;
                if (!auto_y) {
                    ImGui::SetNextItemWidth(100);
                    if (ImGui::InputFloat("Y min", &fy_min, 0, 0, "%.2f")) y_min = fy_min;
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(100);
                    if (ImGui::InputFloat("Y max", &fy_max, 0, 0, "%.2f")) y_max = fy_max;
                }
                if (ImGui::Button("Reset View")) { x_min = -10; x_max = 10; y_min = -10; y_max = 10; }

                ImGui::Separator();
                ImGui::Checkbox("Grid", &show_grid);
                ImGui::SameLine();
                ImGui::Checkbox("Coords", &show_coords);

                ImGui::Separator();
                if (ImGui::Button("Radians")) deg_mode = false;
                ImGui::SameLine();
                if (ImGui::Button("Degrees")) deg_mode = true;
                ImGui::SameLine();
                ImGui::TextDisabled("(%s)", deg_mode ? "deg" : "rad");

                ImGui::EndChild();
                ImGui::SameLine();
            }

            ImVec2 graph_size = ImGui::GetContentRegionAvail();
            if (graph_size.x < 100) graph_size.x = 100;
            if (graph_size.y < 100) graph_size.y = 100;

            ImGui::BeginChild("GraphArea", graph_size, false,
                ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
            ImVec2 canvas_pos = ImGui::GetCursorScreenPos();
            ImVec2 canvas_size = ImGui::GetContentRegionAvail();
            if (canvas_size.x < 10) canvas_size.x = 10;
            if (canvas_size.y < 10) canvas_size.y = 10;

            ImDrawList* dl = ImGui::GetWindowDrawList();
            dl->PushClipRect(canvas_pos, ImVec2(canvas_pos.x + canvas_size.x, canvas_pos.y + canvas_size.y));
            dl->AddRectFilled(canvas_pos, ImVec2(canvas_pos.x + canvas_size.x, canvas_pos.y + canvas_size.y), IM_COL32(20, 20, 28, 255));

            auto w2s = [&](double wx, double wy) -> ImVec2 {
                float sx = canvas_pos.x + (float)((wx - x_min) / (x_max - x_min) * canvas_size.x);
                float sy = canvas_pos.y + (float)((1.0 - (wy - y_min) / (y_max - y_min)) * canvas_size.y);
                return ImVec2(sx, sy);
            };
            auto s2w = [&](float sx, float sy) -> std::pair<double, double> {
                double wx = x_min + (sx - canvas_pos.x) / canvas_size.x * (x_max - x_min);
                double wy = y_max - (sy - canvas_pos.y) / canvas_size.y * (y_max - y_min);
                return {wx, wy};
            };

            if (auto_y) {
                double y_lo = 1e30, y_hi = -1e30;
                bool any_valid = false;
                for (auto& f : funcs) {
                    if (!f.visible || !ensure_ast(f)) continue;
                    if (f.gtype == GT_EXPLICIT && f.cached_ast) {
                        int n = 200;
                        double ddx = (x_max - x_min) / n;
                        for (int i = 0; i <= n; i++) {
                            EvalVars v; v.x = x_min + i * ddx;
                            double yv = fast_eval(f.cached_ast, v, deg_mode);
                            if (std::isfinite(yv) && fabs(yv) < 1e8) {
                                if (yv < y_lo) y_lo = yv;
                                if (yv > y_hi) y_hi = yv;
                                any_valid = true;
                            }
                        }
                    } else if (f.gtype == GT_PARAMETRIC && f.cached_ast && f.cached_ast2) {
                        int n = 500;
                        double dt = (f.t_max - f.t_min) / n;
                        for (int i = 0; i <= n; i++) {
                            EvalVars v; v.t = f.t_min + i * dt;
                            double xv = fast_eval(f.cached_ast, v, deg_mode);
                            double yv = fast_eval(f.cached_ast2, v, deg_mode);
                            if (std::isfinite(xv) && std::isfinite(yv) && fabs(yv) < 1e8) {
                                if (yv < y_lo) y_lo = yv;
                                if (yv > y_hi) y_hi = yv;
                                any_valid = true;
                            }
                        }
                    } else if (f.gtype == GT_POLAR && f.cached_ast) {
                        int n = 500;
                        double dth = (f.t_max - f.t_min) / n;
                        for (int i = 0; i <= n; i++) {
                            EvalVars v; v.theta = f.t_min + i * dth;
                            double r = fast_eval(f.cached_ast, v, deg_mode);
                            if (std::isfinite(r)) {
                                double yv = r * sin(f.t_min + i * dth);
                                if (fabs(yv) < 1e8) {
                                    if (yv < y_lo) y_lo = yv;
                                    if (yv > y_hi) y_hi = yv;
                                    any_valid = true;
                                }
                            }
                        }
                    }
                }
                if (any_valid) {
                    double margin = (y_hi - y_lo) * 0.1;
                    if (margin < 0.5) margin = 0.5;
                    y_min = y_lo - margin;
                    y_max = y_hi + margin;
                }
            }

            {
                double x_range = x_max - x_min;
                double y_range = y_max - y_min;
                double x_ppu = canvas_size.x / x_range;
                double y_ppu = canvas_size.y / y_range;
                if (x_ppu > y_ppu) {
                    double new_x_range = canvas_size.x / y_ppu;
                    double cx = (x_min + x_max) / 2;
                    x_min = cx - new_x_range / 2;
                    x_max = cx + new_x_range / 2;
                } else {
                    double new_y_range = canvas_size.y / x_ppu;
                    double cy = (y_min + y_max) / 2;
                    y_min = cy - new_y_range / 2;
                    y_max = cy + new_y_range / 2;
                }
            }

            if (show_grid) {
                double x_step = nice_step(x_max - x_min);
                double y_step = nice_step(y_max - y_min);
                for (double gx = ceil(x_min / x_step) * x_step; gx <= x_max; gx += x_step) {
                    dl->AddLine(w2s(gx, y_min), w2s(gx, y_max), IM_COL32(50, 50, 60, 255), 1.0f);
                    ImVec2 lp = w2s(gx, y_min);
                    dl->AddText(ImVec2(lp.x + 2, lp.y - 16), IM_COL32(150, 150, 160, 255), format_tick(gx).c_str());
                }
                for (double gy = ceil(y_min / y_step) * y_step; gy <= y_max; gy += y_step) {
                    dl->AddLine(w2s(x_min, gy), w2s(x_max, gy), IM_COL32(50, 50, 60, 255), 1.0f);
                    ImVec2 lp = w2s(x_min, gy);
                    dl->AddText(ImVec2(lp.x + 4, lp.y + 2), IM_COL32(150, 150, 160, 255), format_tick(gy).c_str());
                }
            }

            dl->AddLine(w2s(0, y_min), w2s(0, y_max), IM_COL32(120, 120, 140, 255), 2.0f);
            dl->AddLine(w2s(x_min, 0), w2s(x_max, 0), IM_COL32(120, 120, 140, 255), 2.0f);

            for (auto& f : funcs) {
                if (!f.visible) continue;
                switch (f.gtype) {
                case GT_EXPLICIT:
                    draw_explicit(dl, f, deg_mode, x_min, x_max, y_min, y_max, w2s, canvas_size.x);
                    break;
                case GT_IMPLICIT:
                    draw_implicit(dl, f, deg_mode, x_min, x_max, y_min, y_max, w2s, canvas_size.x, canvas_size.y);
                    break;
                case GT_PARAMETRIC:
                    draw_parametric(dl, f, deg_mode, x_min, x_max, y_min, y_max, w2s);
                    break;
                case GT_POLAR:
                    draw_polar(dl, f, deg_mode, x_min, x_max, y_min, y_max, w2s);
                    break;
                }
            }

            if (show_coords && ImGui::IsWindowHovered()) {
                ImVec2 mouse = ImGui::GetIO().MousePos;
                if (mouse.x >= canvas_pos.x && mouse.x <= canvas_pos.x + canvas_size.x &&
                    mouse.y >= canvas_pos.y && mouse.y <= canvas_pos.y + canvas_size.y) {
                    auto [wx, wy] = s2w(mouse.x, mouse.y);
                    dl->AddLine(ImVec2(mouse.x, canvas_pos.y), ImVec2(mouse.x, canvas_pos.y + canvas_size.y), IM_COL32(200, 200, 200, 80), 1.0f);
                    dl->AddLine(ImVec2(canvas_pos.x, mouse.y), ImVec2(canvas_pos.x + canvas_size.x, mouse.y), IM_COL32(200, 200, 200, 80), 1.0f);
                    char coord_buf[128];
                    snprintf(coord_buf, sizeof(coord_buf), "x: %.4g  y: %.4g", wx, wy);
                    dl->AddText(ImVec2(mouse.x + 12, mouse.y - 20), IM_COL32(255, 255, 200, 255), coord_buf);
                    for (auto& f : funcs) {
                        if (!f.visible || !ensure_ast(f) || f.gtype != GT_EXPLICIT || !f.cached_ast) continue;
                        EvalVars v; v.x = wx;
                        double fv = fast_eval(f.cached_ast, v, deg_mode);
                        if (std::isfinite(fv)) {
                            ImVec2 pt = w2s(wx, fv);
                            if (pt.y >= canvas_pos.y && pt.y <= canvas_pos.y + canvas_size.y) {
                                dl->AddCircleFilled(pt, 5.0f, f.color);
                                char val_buf[128];
                                snprintf(val_buf, sizeof(val_buf), "%.4g", fv);
                                dl->AddText(ImVec2(pt.x + 8, pt.y - 8), f.color, val_buf);
                            }
                        }
                    }
                }
            }

            {
                ImGui::InvisibleButton("graph_canvas", canvas_size);
                if (ImGui::IsItemHovered()) {
                    float wheel = ImGui::GetIO().MouseWheel;
                    if (wheel != 0) {
                        ImVec2 mouse = ImGui::GetIO().MousePos;
                        auto [wx, wy] = s2w(mouse.x, mouse.y);
                        double factor = pow(0.9, wheel);
                        x_min = wx + (x_min - wx) * factor;
                        x_max = wx + (x_max - wx) * factor;
                        y_min = wy + (y_min - wy) * factor;
                        y_max = wy + (y_max - wy) * factor;
                    }
                }
                if (ImGui::IsItemActive() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    dragging = true;
                    drag_start = ImGui::GetIO().MousePos;
                    drag_xmin = x_min; drag_xmax = x_max;
                    drag_ymin = y_min; drag_ymax = y_max;
                }
                if (dragging && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                    ImVec2 delta = ImVec2(ImGui::GetIO().MousePos.x - drag_start.x, ImGui::GetIO().MousePos.y - drag_start.y);
                    double scale = (drag_xmax - drag_xmin) / canvas_size.x;
                    double dxw = -delta.x * scale;
                    double dyw = delta.y * scale;
                    x_min = drag_xmin + dxw; x_max = drag_xmax + dxw;
                    y_min = drag_ymin + dyw; y_max = drag_ymax + dyw;
                }
                if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) dragging = false;
            }

            dl->PopClipRect();
            ImGui::EndChild();

            ImGui::SetNextWindowPos(ImVec2(canvas_pos.x + 8 * main_scale, canvas_pos.y + 8 * main_scale));
            ImGui::SetNextWindowBgAlpha(0.6f);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(6, 4));
            ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 4.0f);
            ImGui::Begin("##SettingsBtn", nullptr,
                ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_AlwaysAutoResize |
                ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
                ImGuiWindowFlags_NoNav);
            if (ImGui::Button(show_settings ? "<< Hide" : "Settings >>")) show_settings = !show_settings;
            ImGui::End();
            ImGui::PopStyleVar(2);

            ImGui::PopStyleVar();
            ImGui::End();
        }

        ImGui::Render();
        const float clear_color_with_alpha[4] = {
            clear_color.x * clear_color.w, clear_color.y * clear_color.w,
            clear_color.z * clear_color.w, clear_color.w
        };
        g_pd3dDeviceContext->OMSetRenderTargets(1, &g_mainRenderTargetView, nullptr);
        g_pd3dDeviceContext->ClearRenderTargetView(g_mainRenderTargetView, clear_color_with_alpha);
        ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());

        HRESULT hr = g_pSwapChain->Present(1, 0);
        g_SwapChainOccluded = (hr == DXGI_STATUS_OCCLUDED);
    }

    ImGui_ImplDX11_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImGui::DestroyContext();

    CleanupDeviceD3D();
    ::DestroyWindow(hwnd);
    ::UnregisterClassW(wc.lpszClassName, wc.hInstance);
    return 0;
}

bool CreateDeviceD3D(HWND hWnd) {
    DXGI_SWAP_CHAIN_DESC sd;
    ZeroMemory(&sd, sizeof(sd));
    sd.BufferCount = 2;
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferDesc.RefreshRate = {60, 1};
    sd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow = hWnd;
    sd.SampleDesc = {1, 0};
    sd.Windowed = TRUE;
    sd.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;

    UINT flags = 0;
    D3D_FEATURE_LEVEL fl;
    const D3D_FEATURE_LEVEL fla[] = {D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_0};
    HRESULT res = D3D11CreateDeviceAndSwapChain(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, flags, fla, 2, D3D11_SDK_VERSION, &sd, &g_pSwapChain, &g_pd3dDevice, &fl, &g_pd3dDeviceContext);
    if (res == DXGI_ERROR_UNSUPPORTED)
        res = D3D11CreateDeviceAndSwapChain(nullptr, D3D_DRIVER_TYPE_WARP, nullptr, flags, fla, 2, D3D11_SDK_VERSION, &sd, &g_pSwapChain, &g_pd3dDevice, &fl, &g_pd3dDeviceContext);
    if (res != S_OK) return false;
    CreateRenderTarget();
    return true;
}

void CleanupDeviceD3D() {
    CleanupRenderTarget();
    if (g_pSwapChain) { g_pSwapChain->Release(); g_pSwapChain = nullptr; }
    if (g_pd3dDeviceContext) { g_pd3dDeviceContext->Release(); g_pd3dDeviceContext = nullptr; }
    if (g_pd3dDevice) { g_pd3dDevice->Release(); g_pd3dDevice = nullptr; }
}

void CreateRenderTarget() {
    ID3D11Texture2D* pBackBuffer;
    g_pSwapChain->GetBuffer(0, IID_PPV_ARGS(&pBackBuffer));
    g_pd3dDevice->CreateRenderTargetView(pBackBuffer, nullptr, &g_mainRenderTargetView);
    pBackBuffer->Release();
}

void CleanupRenderTarget() {
    if (g_mainRenderTargetView) { g_mainRenderTargetView->Release(); g_mainRenderTargetView = nullptr; }
}

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND, UINT, WPARAM, LPARAM);

LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    if (ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam)) return true;
    switch (msg) {
    case WM_SIZE: if (wParam == SIZE_MINIMIZED) return 0; g_ResizeWidth = (UINT)LOWORD(lParam); g_ResizeHeight = (UINT)HIWORD(lParam); return 0;
    case WM_SYSCOMMAND: if ((wParam & 0xfff0) == SC_KEYMENU) return 0; break;
    case WM_DESTROY: ::PostQuitMessage(0); return 0;
    }
    return ::DefWindowProcW(hWnd, msg, wParam, lParam);
}
