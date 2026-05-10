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

struct FuncEntry {
    char expr[256];
    bool visible;
    ImU32 color;
    bool valid;
    std::string error_msg;
    ASTPtr cached_ast;
    std::string cached_expr_str;
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
    std::string cur(f.expr);
    if (cur == f.cached_expr_str && f.cached_ast)
        return f.valid;
    f.cached_expr_str = cur;
    f.cached_ast.reset();
    if (cur.empty()) {
        f.valid = true;
        f.error_msg = "";
        return true;
    }
    try {
        Lexer lexer(cur);
        auto tokens = lexer.tokenize();
        Parser parser(tokens);
        f.cached_ast = parser.parse();
        f.valid = true;
        f.error_msg = "";
    } catch (const std::exception& e) {
        f.valid = false;
        f.error_msg = e.what();
    }
    return f.valid;
}

static double fast_eval(ASTPtr node, double x, bool deg_mode) {
    if (!node) return NAN;
    switch (node->type) {
    case ASTNode::NUMBER: {
        double v = std::stod(node->number.to_string());
        return v;
    }
    case ASTNode::CONSTANT:
        if (node->name == "pi") return M_PI;
        if (node->name == "e") return M_E;
        return NAN;
    case ASTNode::VARIABLE:
        if (node->name == "x") return x;
        if (node->name == "ans") return 0;
        return 0;
    case ASTNode::BINOP: {
        double l = fast_eval(node->left, x, deg_mode);
        double r = fast_eval(node->right, x, deg_mode);
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
        if (node->op == '-') return -fast_eval(node->left, x, deg_mode);
        return fast_eval(node->left, x, deg_mode);
    case ASTNode::FACTORIAL: {
        double v = fast_eval(node->left, x, deg_mode);
        if (v < 0 || v != floor(v) || v > 170) return NAN;
        double r = 1;
        for (int i = 2; i <= (int)v; i++) r *= i;
        return r;
    }
    case ASTNode::FUNCTION: {
        auto& args = node->args;
        auto fe = [&](int i) { return fast_eval(args[i], x, deg_mode); };
        const std::string& name = node->name;
        auto to_rad = [&](double v) { return deg_mode ? v * M_PI / 180.0 : v; };
        auto from_rad = [&](double v) { return deg_mode ? v * 180.0 / M_PI : v; };

        if (name == "sin" && args.size() == 1) return sin(to_rad(fe(0)));
        if (name == "cos" && args.size() == 1) return cos(to_rad(fe(0)));
        if (name == "tan" && args.size() == 1) { double v = cos(to_rad(fe(0))); return fabs(v) < 1e-15 ? NAN : tan(to_rad(fe(0))); }
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
        if (name == "log" && args.size() == 2) { double b = fe(0); double v = fe(1); return (b <= 0 || v <= 0) ? NAN : log(v) / log(b); }
        if (name == "exp" && args.size() == 1) return exp(fe(0));
        if (name == "fact" && args.size() == 1) { double v = fe(0); if (v < 0 || v != floor(v) || v > 170) return NAN; double r = 1; for (int i = 2; i <= (int)v; i++) r *= i; return r; }
        if (name == "floor" && args.size() == 1) return floor(fe(0));
        if (name == "ceil" && args.size() == 1) return ceil(fe(0));
        if (name == "round" && args.size() == 1) return round(fe(0));
        if (name == "sign" && args.size() == 1) { double v = fe(0); return v > 0 ? 1 : (v < 0 ? -1 : 0); }
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

int main(int argc, char** argv) {
    std::string initial_expr;
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            if (i > 1) initial_expr += " ";
            initial_expr += argv[i];
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
    if (!initial_expr.empty()) {
        int len = MultiByteToWideChar(CP_UTF8, 0, initial_expr.c_str(), -1, nullptr, 0);
        std::wstring wexpr(len, 0);
        MultiByteToWideChar(CP_UTF8, 0, initial_expr.c_str(), -1, &wexpr[0], len);
        wtitle = L"y = " + wexpr;
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
    if (!initial_expr.empty()) {
        FuncEntry fe;
        memset(fe.expr, 0, sizeof(fe.expr));
        strncpy(fe.expr, initial_expr.c_str(), sizeof(fe.expr) - 1);
        fe.visible = true;
        fe.color = func_colors[0];
        fe.valid = true;
        fe.error_msg = "";
        fe.cached_ast = nullptr;
        fe.cached_expr_str = "";
        funcs.push_back(fe);
    } else {
        funcs.push_back({"sin(x)", true, func_colors[0], true, "", nullptr, ""});
    }
    for (int i = 1; i < 8; i++) {
        funcs.push_back({"", true, func_colors[i % 8], true, "", nullptr, ""});
    }

    double x_min = -10.0, x_max = 10.0;
    double y_min = -5.0, y_max = 5.0;
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
            if (msg.message == WM_QUIT)
                done = true;
        }
        if (done) break;

        if (g_SwapChainOccluded &&
            g_pSwapChain->Present(0, DXGI_PRESENT_TEST) == DXGI_STATUS_OCCLUDED) {
            ::Sleep(10);
            continue;
        }
        g_SwapChainOccluded = false;

        if (g_ResizeWidth != 0 && g_ResizeHeight != 0) {
            CleanupRenderTarget();
            g_pSwapChain->ResizeBuffers(0, g_ResizeWidth, g_ResizeHeight,
                DXGI_FORMAT_UNKNOWN, 0);
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
                float panel_width = 300.0f * main_scale;

                ImGui::BeginChild("SettingsPanel", ImVec2(panel_width, 0), true);

                ImGui::Text("Functions:");
                for (int i = 0; i < (int)funcs.size(); i++) {
                    ImGui::PushID(i);
                    ImGui::Checkbox("##vis", &funcs[i].visible);
                    ImGui::SameLine();
                    ImGui::ColorButton("##col",
                        ImVec4(
                            ((funcs[i].color >> 0) & 0xFF) / 255.0f,
                            ((funcs[i].color >> 8) & 0xFF) / 255.0f,
                            ((funcs[i].color >> 16) & 0xFF) / 255.0f,
                            1.0f),
                        ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoPicker,
                        ImVec2(16, 16));
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(-40);
                    ImGui::InputText("##expr", funcs[i].expr,
                        sizeof(funcs[i].expr));
                    ImGui::SameLine();
                    if (ImGui::Button("X##del")) {
                        funcs.erase(funcs.begin() + i);
                        i--;
                        ImGui::PopID();
                        continue;
                    }
                    ensure_ast(funcs[i]);
                    if (!funcs[i].valid && funcs[i].expr[0] != '\0') {
                        ImGui::TextColored(ImVec4(1, 0.5f, 0.5f, 1), "  %s",
                            funcs[i].error_msg.c_str());
                    }
                    ImGui::PopID();
                }
                if (ImGui::Button("+ Add")) {
                    int ci = (int)funcs.size() % 8;
                    funcs.push_back({"", true, func_colors[ci], true, "", nullptr, ""});
                }

                ImGui::Separator();
                ImGui::Checkbox("Auto Y range", &auto_y);
                float fx_min = (float)x_min, fx_max = (float)x_max;
                float fy_min = (float)y_min, fy_max = (float)y_max;
                ImGui::SetNextItemWidth(100);
                if (ImGui::InputFloat("X min", &fx_min, 0, 0, "%.2f"))
                    x_min = fx_min;
                ImGui::SameLine();
                ImGui::SetNextItemWidth(100);
                if (ImGui::InputFloat("X max", &fx_max, 0, 0, "%.2f"))
                    x_max = fx_max;
                if (!auto_y) {
                    ImGui::SetNextItemWidth(100);
                    if (ImGui::InputFloat("Y min", &fy_min, 0, 0, "%.2f"))
                        y_min = fy_min;
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(100);
                    if (ImGui::InputFloat("Y max", &fy_max, 0, 0, "%.2f"))
                        y_max = fy_max;
                }
                if (ImGui::Button("Reset View")) {
                    x_min = -10; x_max = 10;
                    y_min = -5; y_max = 5;
                }

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
            dl->PushClipRect(canvas_pos,
                ImVec2(canvas_pos.x + canvas_size.x,
                       canvas_pos.y + canvas_size.y));

            dl->AddRectFilled(canvas_pos,
                ImVec2(canvas_pos.x + canvas_size.x,
                       canvas_pos.y + canvas_size.y),
                IM_COL32(20, 20, 28, 255));

            auto world_to_screen = [&](double wx, double wy) -> ImVec2 {
                float sx = canvas_pos.x + (float)((wx - x_min) / (x_max - x_min) * canvas_size.x);
                float sy = canvas_pos.y + (float)((1.0 - (wy - y_min) / (y_max - y_min)) * canvas_size.y);
                return ImVec2(sx, sy);
            };

            auto screen_to_world = [&](float sx, float sy) -> std::pair<double, double> {
                double wx = x_min + (sx - canvas_pos.x) / canvas_size.x * (x_max - x_min);
                double wy = y_max - (sy - canvas_pos.y) / canvas_size.y * (y_max - y_min);
                return {wx, wy};
            };

            if (auto_y) {
                double y_lo = 1e30, y_hi = -1e30;
                bool any_valid = false;
                int n_auto = 200;
                double dx_auto = (x_max - x_min) / n_auto;
                for (auto& f : funcs) {
                    if (!f.visible || !ensure_ast(f) || !f.cached_ast) continue;
                    for (int i = 0; i <= n_auto; i++) {
                        double xv = x_min + i * dx_auto;
                        double yv = fast_eval(f.cached_ast, xv, deg_mode);
                        if (std::isfinite(yv) && fabs(yv) < 1e8) {
                            if (yv < y_lo) y_lo = yv;
                            if (yv > y_hi) y_hi = yv;
                            any_valid = true;
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

            if (show_grid) {
                double x_step = nice_step(x_max - x_min);
                double y_step = nice_step(y_max - y_min);

                double x_start = ceil(x_min / x_step) * x_step;
                for (double gx = x_start; gx <= x_max; gx += x_step) {
                    ImVec2 p1 = world_to_screen(gx, y_min);
                    ImVec2 p2 = world_to_screen(gx, y_max);
                    dl->AddLine(p1, p2, IM_COL32(50, 50, 60, 255), 1.0f);
                    ImVec2 label_pos = world_to_screen(gx, y_min);
                    dl->AddText(ImVec2(label_pos.x + 2, label_pos.y - 16),
                        IM_COL32(150, 150, 160, 255), format_tick(gx).c_str());
                }

                double y_start = ceil(y_min / y_step) * y_step;
                for (double gy = y_start; gy <= y_max; gy += y_step) {
                    ImVec2 p1 = world_to_screen(x_min, gy);
                    ImVec2 p2 = world_to_screen(x_max, gy);
                    dl->AddLine(p1, p2, IM_COL32(50, 50, 60, 255), 1.0f);
                    ImVec2 label_pos = world_to_screen(x_min, gy);
                    dl->AddText(ImVec2(label_pos.x + 4, label_pos.y + 2),
                        IM_COL32(150, 150, 160, 255), format_tick(gy).c_str());
                }
            }

            {
                ImVec2 p1 = world_to_screen(0, y_min);
                ImVec2 p2 = world_to_screen(0, y_max);
                dl->AddLine(p1, p2, IM_COL32(120, 120, 140, 255), 2.0f);
                ImVec2 p3 = world_to_screen(x_min, 0);
                ImVec2 p4 = world_to_screen(x_max, 0);
                dl->AddLine(p3, p4, IM_COL32(120, 120, 140, 255), 2.0f);
            }

            for (auto& f : funcs) {
                if (!f.visible || !ensure_ast(f) || !f.cached_ast) continue;

                int n_samples = (int)canvas_size.x;
                if (n_samples < 100) n_samples = 100;
                if (n_samples > 2000) n_samples = 2000;
                double dx = (x_max - x_min) / n_samples;

                std::vector<ImVec2> points;
                points.reserve(n_samples + 1);
                double prev_y = NAN;

                for (int i = 0; i <= n_samples; i++) {
                    double xv = x_min + i * dx;
                    double yv = fast_eval(f.cached_ast, xv, deg_mode);

                    if (!std::isfinite(yv) || fabs(yv) > 1e10) {
                        if (points.size() >= 2) {
                            dl->AddPolyline(points.data(), (int)points.size(),
                                f.color, false, 2.0f);
                        }
                        points.clear();
                        prev_y = NAN;
                        continue;
                    }

                    if (std::isfinite(prev_y) && fabs(yv - prev_y) > (y_max - y_min) * 10) {
                        if (points.size() >= 2) {
                            dl->AddPolyline(points.data(), (int)points.size(),
                                f.color, false, 2.0f);
                        }
                        points.clear();
                    }

                    ImVec2 sp = world_to_screen(xv, yv);
                    points.push_back(sp);
                    prev_y = yv;
                }
                if (points.size() >= 2) {
                    dl->AddPolyline(points.data(), (int)points.size(),
                        f.color, false, 2.0f);
                }
            }

            if (show_coords && ImGui::IsWindowHovered()) {
                ImVec2 mouse = ImGui::GetIO().MousePos;
                if (mouse.x >= canvas_pos.x && mouse.x <= canvas_pos.x + canvas_size.x &&
                    mouse.y >= canvas_pos.y && mouse.y <= canvas_pos.y + canvas_size.y) {
                    auto [wx, wy] = screen_to_world(mouse.x, mouse.y);

                    dl->AddLine(ImVec2(mouse.x, canvas_pos.y),
                        ImVec2(mouse.x, canvas_pos.y + canvas_size.y),
                        IM_COL32(200, 200, 200, 80), 1.0f);
                    dl->AddLine(ImVec2(canvas_pos.x, mouse.y),
                        ImVec2(canvas_pos.x + canvas_size.x, mouse.y),
                        IM_COL32(200, 200, 200, 80), 1.0f);

                    char coord_buf[128];
                    snprintf(coord_buf, sizeof(coord_buf), "x: %.4g  y: %.4g", wx, wy);
                    dl->AddText(ImVec2(mouse.x + 12, mouse.y - 20),
                        IM_COL32(255, 255, 200, 255), coord_buf);

                    for (auto& f : funcs) {
                        if (!f.visible || !ensure_ast(f) || !f.cached_ast) continue;
                        double fv = fast_eval(f.cached_ast, wx, deg_mode);
                        if (std::isfinite(fv)) {
                            ImVec2 pt = world_to_screen(wx, fv);
                            if (pt.y >= canvas_pos.y && pt.y <= canvas_pos.y + canvas_size.y) {
                                dl->AddCircleFilled(pt, 5.0f, f.color);
                                char val_buf[128];
                                snprintf(val_buf, sizeof(val_buf), "%.4g", fv);
                                dl->AddText(ImVec2(pt.x + 8, pt.y - 8),
                                    f.color, val_buf);
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
                        auto [wx, unused] = screen_to_world(mouse.x, mouse.y);
                        (void)unused;
                        double factor = pow(0.9, wheel);
                        x_min = wx + (x_min - wx) * factor;
                        x_max = wx + (x_max - wx) * factor;
                    }
                }
                if (ImGui::IsItemActive() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    dragging = true;
                    drag_start = ImGui::GetIO().MousePos;
                    drag_xmin = x_min; drag_xmax = x_max;
                    drag_ymin = y_min; drag_ymax = y_max;
                }
                if (dragging && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                    ImVec2 delta = ImVec2(
                        ImGui::GetIO().MousePos.x - drag_start.x,
                        ImGui::GetIO().MousePos.y - drag_start.y);
                    double dx_world = -delta.x / canvas_size.x * (drag_xmax - drag_xmin);
                    double dy_world = delta.y / canvas_size.y * (drag_ymax - drag_ymin);
                    x_min = drag_xmin + dx_world;
                    x_max = drag_xmax + dx_world;
                    y_min = drag_ymin + dy_world;
                    y_max = drag_ymax + dy_world;
                }
                if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                    dragging = false;
                }
            }

            dl->PopClipRect();
            ImGui::EndChild();

            ImGui::SetNextWindowPos(ImVec2(canvas_pos.x + 8 * main_scale,
                                            canvas_pos.y + 8 * main_scale));
            ImGui::SetNextWindowBgAlpha(0.6f);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(6, 4));
            ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 4.0f);
            ImGui::Begin("##SettingsBtn", nullptr,
                ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_AlwaysAutoResize |
                ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
                ImGuiWindowFlags_NoNav);
            if (ImGui::Button(show_settings ? "<< Hide" : "Settings >>")) {
                show_settings = !show_settings;
            }
            ImGui::End();
            ImGui::PopStyleVar(2);

            ImGui::PopStyleVar();
            ImGui::End();
        }

        ImGui::Render();
        const float clear_color_with_alpha[4] = {
            clear_color.x * clear_color.w,
            clear_color.y * clear_color.w,
            clear_color.z * clear_color.w,
            clear_color.w
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
    sd.BufferDesc.Width = 0;
    sd.BufferDesc.Height = 0;
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferDesc.RefreshRate.Numerator = 60;
    sd.BufferDesc.RefreshRate.Denominator = 1;
    sd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow = hWnd;
    sd.SampleDesc.Count = 1;
    sd.SampleDesc.Quality = 0;
    sd.Windowed = TRUE;
    sd.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;

    UINT createDeviceFlags = 0;
    D3D_FEATURE_LEVEL featureLevel;
    const D3D_FEATURE_LEVEL featureLevelArray[2] = {
        D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_0 };
    HRESULT res = D3D11CreateDeviceAndSwapChain(nullptr,
        D3D_DRIVER_TYPE_HARDWARE, nullptr, createDeviceFlags,
        featureLevelArray, 2, D3D11_SDK_VERSION, &sd,
        &g_pSwapChain, &g_pd3dDevice, &featureLevel,
        &g_pd3dDeviceContext);
    if (res == DXGI_ERROR_UNSUPPORTED)
        res = D3D11CreateDeviceAndSwapChain(nullptr,
            D3D_DRIVER_TYPE_WARP, nullptr, createDeviceFlags,
            featureLevelArray, 2, D3D11_SDK_VERSION, &sd,
            &g_pSwapChain, &g_pd3dDevice, &featureLevel,
            &g_pd3dDeviceContext);
    if (res != S_OK)
        return false;

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
    if (g_mainRenderTargetView) {
        g_mainRenderTargetView->Release();
        g_mainRenderTargetView = nullptr;
    }
}

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(
    HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    if (ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam))
        return true;

    switch (msg) {
    case WM_SIZE:
        if (wParam == SIZE_MINIMIZED)
            return 0;
        g_ResizeWidth = (UINT)LOWORD(lParam);
        g_ResizeHeight = (UINT)HIWORD(lParam);
        return 0;
    case WM_SYSCOMMAND:
        if ((wParam & 0xfff0) == SC_KEYMENU)
            return 0;
        break;
    case WM_DESTROY:
        ::PostQuitMessage(0);
        return 0;
    }
    return ::DefWindowProcW(hWnd, msg, wParam, lParam);
}
