#include "ccalc.h"
#include <fstream>
#include <conio.h>
#ifdef _WIN32
#include <windows.h>
#endif

static bool is_assignable_var(const std::string& name) {
    return name == "x" || name == "y" ||
           name == "A" || name == "B" || name == "C" || name == "D" ||
           name == "P" ||
           name == "VerA" || name == "VerB" || name == "VerC" || name == "VerD";
}

static bool is_vec_var(const std::string& name) {
    return name == "VerA" || name == "VerB" || name == "VerC" || name == "VerD";
}

static const std::vector<std::string> all_function_names = {
    "sin", "cos", "tan", "asin", "acos", "atan",
    "sinh", "cosh", "tanh",
    "sqrt", "cbrt", "nrt",
    "abs", "ln", "lg", "log",
    "fact", "perm", "P", "comb", "C",
    "gcd", "lcm", "factor",
    "int", "diff", "sum", "prod",
    "complex", "re", "im", "conj", "arg", "mod",
    "rand", "randint",
    "convert", "solve",
    "deg", "rad",
    "vecmod", "dot", "cross", "scalarmul", "mixed", "proj", "decompose",
    "matrix", "det", "inv", "eigen", "trace", "transpose", "identity",
    "mean", "stddev", "variance", "median",
    "simplify", "exp", "floor", "ceil", "round", "sign", "max", "min",
    "pow", "root", "log2", "atan2", "hypot",
    "taylor", "limit", "inttable", "recur"
};

static bool is_known_func(const std::string& name) {
    for (auto& f : all_function_names) {
        if (f == name) return true;
    }
    return false;
}

static void print_help() {
    std::cout << "CCalc - Command-line Advanced Calculator" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Arithmetic:  +  -  *  /  ^  %" << std::endl;
    std::cout << "Constants:   pi, e" << std::endl;
    std::cout << "Variables:   x, y, A, B, C, D (default 0, assign: x=expr)" << std::endl;
    std::cout << "Vectors:     VerA, VerB, VerC, VerD (default empty, assign: VerA=(1,2,3))" << std::endl;
    std::cout << "Base input:  0b1010(binary) 0o77(octal) 0xFF(hex)" << std::endl;
    std::cout << std::endl;
    std::cout << "Functions:" << std::endl;
    std::cout << "  sin(x), cos(x), tan(x)          - Trigonometric" << std::endl;
    std::cout << "  asin(x), acos(x), atan(x)       - Inverse trig" << std::endl;
    std::cout << "  sinh(x), cosh(x), tanh(x)       - Hyperbolic" << std::endl;
    std::cout << "  sqrt(x), cbrt(x), nrt(n, x)     - Roots" << std::endl;
    std::cout << "  abs(x)                          - Absolute value" << std::endl;
    std::cout << "  ln(x), lg(x), log(base, x)      - Logarithms" << std::endl;
    std::cout << "  exp(x)                          - Exponential" << std::endl;
    std::cout << "  n! or fact(n)                   - Factorial" << std::endl;
    std::cout << "  P(n,k) or perm(n,k)             - Permutation" << std::endl;
    std::cout << "  C(n,k) or comb(n,k)             - Combination" << std::endl;
    std::cout << "  gcd(a,b), lcm(a,b)              - GCD/LCM" << std::endl;
    std::cout << "  factor(n)                       - Prime factorization" << std::endl;
    std::cout << "  floor(x), ceil(x), round(x)     - Rounding" << std::endl;
    std::cout << "  sign(x), max(a,b), min(a,b)     - Other" << std::endl;
    std::cout << std::endl;
    std::cout << "Calculus:" << std::endl;
    std::cout << "  int(expr, var, a, b)            - Definite integral" << std::endl;
    std::cout << "  diff(expr, var, point)          - Derivative at point" << std::endl;
    std::cout << "  diff(expr, var)                 - Symbolic derivative" << std::endl;
    std::cout << "  taylor(expr, var, point, n)     - Taylor expansion" << std::endl;
    std::cout << "  limit(expr, var, val)           - Numerical limit" << std::endl;
    std::cout << "  inttable(expr)                  - Integration table lookup" << std::endl;
    std::cout << "  sum(expr, var, from, to)        - Summation" << std::endl;
    std::cout << "  prod(expr, var, from, to)       - Product" << std::endl;
    std::cout << std::endl;
    std::cout << "Sequences:" << std::endl;
    std::cout << "  recur(expr, var, init, n)       - Recursive sequence" << std::endl;
    std::cout << std::endl;
    std::cout << "Custom functions:" << std::endl;
    std::cout << "  f(x) := x^2+1                  - Define custom function" << std::endl;
    std::cout << std::endl;
    std::cout << "Complex:" << std::endl;
    std::cout << "  complex(a, b)                   - Complex number a+bi" << std::endl;
    std::cout << "  re(z), im(z), conj(z)           - Complex operations" << std::endl;
    std::cout << "  arg(z), mod(z)                  - Argument, modulus" << std::endl;
    std::cout << std::endl;
    std::cout << "Other:" << std::endl;
    std::cout << "  rand()                          - Random [0,1)" << std::endl;
    std::cout << "  randint(a, b)                   - Random integer" << std::endl;
    std::cout << "  convert(val, from, to)          - Unit conversion" << std::endl;
    std::cout << std::endl;
    std::cout << "Vector:" << std::endl;
    std::cout << "  vecmod(v)                       - Vector magnitude" << std::endl;
    std::cout << "  dot(a, b)                       - Dot product" << std::endl;
    std::cout << "  cross(a, b)                     - Cross product (3D)" << std::endl;
    std::cout << "  scalarmul(s, v)                 - Scalar multiplication" << std::endl;
    std::cout << "  mixed(a, b, c)                  - Mixed/scalar triple product (3D)" << std::endl;
    std::cout << "  proj(a, b)                      - Projection of a onto b" << std::endl;
    std::cout << "  decompose(a, b, c)              - Decompose 2D: a=alpha*b+beta*c" << std::endl;
    std::cout << "  decompose(a, b, c, d)           - Decompose 3D: a=alpha*b+beta*c+gamma*d" << std::endl;
    std::cout << "Equation:" << std::endl;
    std::cout << "  solve(a,b,c)                    - Quadratic: ax^2+bx+c=0 (exact)" << std::endl;
    std::cout << "  solve(a,b,c,d)                  - Cubic: ax^3+bx^2+cx+d=0" << std::endl;
    std::cout << "  solve(a,b,c,d,e)                - Quartic: ax^4+bx^3+cx^2+dx+e=0" << std::endl;
    std::cout << "  graph(expr)                     - Plot y=expr (e.g. graph(sin(x)))" << std::endl;
    std::cout << "  graph_implicit(expr)            - Plot f(x,y)=0 (e.g. graph_implicit(x^2+y^2-1))" << std::endl;
    std::cout << "  graph_param(t, x_expr, y_expr)  - Plot parametric (e.g. graph_param(t, cos(t), sin(t)))" << std::endl;
    std::cout << "  graph_polar(r_expr)             - Plot polar r=f(theta) (e.g. graph_polar(1+cos(theta)))" << std::endl;
    std::cout << std::endl;
    std::cout << "Matrix:" << std::endl;
    std::cout << "  [[1,2],[3,4]]                  - Matrix literal" << std::endl;
    std::cout << "  det(m), inv(m), eigen(m)       - Determinant, inverse, eigenvalues" << std::endl;
    std::cout << "  trace(m), transpose(m)         - Trace, transpose" << std::endl;
    std::cout << "  identity(n)                    - n x n identity matrix" << std::endl;
    std::cout << std::endl;
    std::cout << "Statistics:" << std::endl;
    std::cout << "  mean(a,b,...)                  - Arithmetic mean" << std::endl;
    std::cout << "  stddev(a,b,...)                - Standard deviation" << std::endl;
    std::cout << "  variance(a,b,...)              - Variance" << std::endl;
    std::cout << "  median(a,b,...)                - Median" << std::endl;
    std::cout << std::endl;
    std::cout << "  simplify(expr)                 - Simplify expression" << std::endl;
    std::cout << std::endl;
    std::cout << "Commands:" << std::endl;
    std::cout << "  deg / rad                       - Set angle mode" << std::endl;
    std::cout << "  base N                          - Set output base (2-36)" << std::endl;
    std::cout << "  prec N                          - Set precision (digits)" << std::endl;
    std::cout << "  latex                           - Toggle LaTeX output mode" << std::endl;
    std::cout << "  ans                             - Last answer" << std::endl;
    std::cout << "  help                            - Show this help" << std::endl;
    std::cout << "  quit / exit                     - Exit program" << std::endl;
    std::cout << std::endl;
    std::cout << "REPL:" << std::endl;
    std::cout << "  Up/Down                         - Navigate history" << std::endl;
    std::cout << "  Tab                             - Complete function name" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  sin(pi/2)        => 1" << std::endl;
    std::cout << "  cos(pi/3)        => 1/2" << std::endl;
    std::cout << "  diff(x^2, x)     => 2*x" << std::endl;
    std::cout << "  taylor(sin(x), x, 0, 5)  => x - x^3/6 + x^5/120" << std::endl;
    std::cout << "  f(x) := x^2+1    => defines f" << std::endl;
    std::cout << "  f(3)             => 10" << std::endl;
    std::cout << "  0b101+0b110      => 11" << std::endl;
}

static std::string trim(const std::string& s) {
    size_t start = 0;
    while (start < s.size()) {
        unsigned char c = s[start];
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') { start++; continue; }
        if (c == 0xEF && start + 2 < s.size() &&
            (unsigned char)s[start+1] == 0xBB && (unsigned char)s[start+2] == 0xBF) {
            start += 3; continue;
        }
        break;
    }
    size_t end = s.size();
    while (end > start) {
        unsigned char c = s[end - 1];
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') { end--; continue; }
        break;
    }
    return s.substr(start, end - start);
}

static std::vector<std::string> split_top_level(const std::string& input) {
    std::vector<std::string> parts;
    int depth = 0;
    int bracket_depth = 0;
    size_t start = 0;
    for (size_t i = 0; i < input.size(); i++) {
        char c = input[i];
        if (c == '(') depth++;
        else if (c == ')') depth--;
        else if (c == '[') bracket_depth++;
        else if (c == ']') bracket_depth--;
        else if (c == ',' && depth == 0 && bracket_depth == 0) {
            parts.push_back(trim(input.substr(start, i - start)));
            start = i + 1;
        }
    }
    parts.push_back(trim(input.substr(start)));
    return parts;
}

struct GraphCmd {
    enum Type { EXPLICIT, IMPLICIT, PARAMETRIC, POLAR, NONE };
    Type type = NONE;
    std::string expr, expr2;
};

static std::vector<std::string> split_args(const std::string& input) {
    std::vector<std::string> parts;
    int depth = 0;
    size_t start = 0;
    for (size_t i = 0; i < input.size(); i++) {
        char c = input[i];
        if (c == '(') depth++;
        else if (c == ')') depth--;
        else if (c == ',' && depth == 0) {
            parts.push_back(trim(input.substr(start, i - start)));
            start = i + 1;
        }
    }
    parts.push_back(trim(input.substr(start)));
    return parts;
}

static GraphCmd parse_graph_cmd(const std::string& input) {
    std::string s = trim(input);
    GraphCmd cmd;

    if (s.size() >= 7 && s.substr(0, 6) == "graph(" && s.back() == ')') {
        cmd.type = GraphCmd::EXPLICIT;
        cmd.expr = trim(s.substr(6, s.size() - 7));
    } else if (s.size() >= 17 && s.substr(0, 15) == "graph_implicit(" && s.back() == ')') {
        cmd.type = GraphCmd::IMPLICIT;
        cmd.expr = trim(s.substr(15, s.size() - 16));
    } else if (s.size() >= 13 && s.substr(0, 11) == "graph_param(" && s.back() == ')') {
        cmd.type = GraphCmd::PARAMETRIC;
        std::string inner = s.substr(11, s.size() - 12);
        auto parts = split_args(inner);
        if (parts.size() >= 3) {
            cmd.expr = trim(parts[1]);
            cmd.expr2 = trim(parts[2]);
        }
    } else if (s.size() >= 13 && s.substr(0, 11) == "graph_polar(" && s.back() == ')') {
        cmd.type = GraphCmd::POLAR;
        cmd.expr = trim(s.substr(11, s.size() - 12));
    }

    return cmd;
}

static std::string syntax_highlight(const std::string& input) {
    std::string result;
    size_t i = 0;
    while (i < input.size()) {
        if (std::isdigit(input[i]) || (input[i] == '.' && i + 1 < input.size() && std::isdigit(input[i+1]))) {
            if (i + 1 < input.size() && (input[i+1] == 'b' || input[i+1] == 'B' ||
                input[i+1] == 'o' || input[i+1] == 'O' ||
                input[i+1] == 'x' || input[i+1] == 'X') && input[i] == '0') {
                std::string num;
                num += input[i]; i++;
                num += input[i]; i++;
                while (i < input.size() && (std::isxdigit(input[i]) || input[i] == '.')) {
                    num += input[i]; i++;
                }
                result += "\033[90m" + num + "\033[0m";
            } else {
                std::string num;
                while (i < input.size() && (std::isdigit(input[i]) || input[i] == '.')) {
                    num += input[i]; i++;
                }
                result += "\033[90m" + num + "\033[0m";
            }
        } else if (std::isalpha(input[i]) || input[i] == '_') {
            std::string ident;
            while (i < input.size() && (std::isalnum(input[i]) || input[i] == '_')) {
                ident += input[i]; i++;
            }
            if (is_known_func(ident)) {
                result += "\033[32m" + ident + "\033[0m";
            } else if (ident == "pi" || ident == "PI" || ident == "e" || ident == "E") {
                result += "\033[90m" + ident + "\033[0m";
            } else {
                result += ident;
            }
        } else if (input[i] == '+' || input[i] == '-' || input[i] == '*' ||
                   input[i] == '/' || input[i] == '^' || input[i] == '%') {
            result += "\033[33m" + std::string(1, input[i]) + "\033[0m";
            i++;
        } else {
            result += input[i];
            i++;
        }
    }
    return result;
}

static std::string readline_custom(const std::string& prompt,
                                    std::vector<std::string>& history,
                                    int& history_idx) {
    std::string line;
    int cursor = 0;

    std::cout << prompt << std::flush;

    while (true) {
        int ch = _getch();

        if (ch == '\r' || ch == '\n') {
            std::cout << std::endl;
            if (!line.empty()) {
                bool found = false;
                for (auto& h : history) {
                    if (h == line) { found = true; break; }
                }
                if (!found) history.push_back(line);
            }
            history_idx = (int)history.size();
            return line;
        }

        if (ch == '\b' || ch == 127) {
            if (cursor > 0) {
                line.erase(cursor - 1, 1);
                cursor--;
                std::cout << "\r" << prompt << syntax_highlight(line)
                          << std::string(line.size() - cursor + 1, ' ');
                if (line.size() - cursor + 1 > 0) {
                    std::cout << "\033[" << (line.size() - cursor + 1) << "D";
                }
                std::cout << std::flush;
                int target = cursor + (int)prompt.size();
                std::cout << "\r\033[" << target << "C" << std::flush;
            }
            continue;
        }

        if (ch == 0x00 || ch == 0xE0) {
            int ch2 = _getch();
            if (ch2 == 0x48) {
                if (history_idx > 0) {
                    history_idx--;
                    line = history[history_idx];
                    cursor = (int)line.size();
                    std::cout << "\r" << prompt << syntax_highlight(line)
                              << std::string(20, ' ') << "\r";
                    std::cout << prompt << syntax_highlight(line) << std::flush;
                }
            } else if (ch2 == 0x50) {
                if (history_idx < (int)history.size() - 1) {
                    history_idx++;
                    line = history[history_idx];
                } else {
                    history_idx = (int)history.size();
                    line.clear();
                }
                cursor = (int)line.size();
                std::cout << "\r" << prompt << syntax_highlight(line)
                          << std::string(20, ' ') << "\r";
                std::cout << prompt << syntax_highlight(line) << std::flush;
            } else if (ch2 == 0x4D) {
                if (cursor < (int)line.size()) {
                    cursor++;
                    std::cout << "\033[C" << std::flush;
                }
            } else if (ch2 == 0x4B) {
                if (cursor > 0) {
                    cursor--;
                    std::cout << "\033[D" << std::flush;
                }
            }
            continue;
        }

        if (ch == '\t') {
            std::string word;
            int ws = cursor - 1;
            while (ws >= 0 && (std::isalnum(line[ws]) || line[ws] == '_')) ws--;
            ws++;
            word = line.substr(ws, cursor - ws);
            if (!word.empty()) {
                std::vector<std::string> matches;
                for (auto& fn : all_function_names) {
                    if (fn.size() >= word.size() && fn.substr(0, word.size()) == word) {
                        matches.push_back(fn);
                    }
                }
                if (matches.size() == 1) {
                    std::string completion = matches[0].substr(word.size());
                    line.insert(cursor, completion);
                    cursor += (int)completion.size();
                    std::cout << "\r" << prompt << syntax_highlight(line)
                              << std::string(20, ' ') << "\r";
                    std::cout << prompt << syntax_highlight(line) << std::flush;
                } else if (matches.size() > 1) {
                    std::string common_prefix = matches[0];
                    for (auto& m : matches) {
                        size_t j = 0;
                        while (j < common_prefix.size() && j < m.size() && common_prefix[j] == m[j]) j++;
                        common_prefix = common_prefix.substr(0, j);
                    }
                    if (common_prefix.size() > word.size()) {
                        std::string completion = common_prefix.substr(word.size());
                        line.insert(cursor, completion);
                        cursor += (int)completion.size();
                    }
                    std::cout << "\r" << prompt << syntax_highlight(line)
                              << std::string(20, ' ') << "\r";
                    std::cout << prompt << syntax_highlight(line) << std::flush;
                }
            }
            continue;
        }

        if (ch >= 32 && ch < 127) {
            line.insert(cursor, 1, (char)ch);
            cursor++;
            std::cout << "\r" << prompt << syntax_highlight(line)
                      << std::string(20, ' ') << "\r";
            std::cout << prompt << syntax_highlight(line) << std::flush;
            int target = cursor + (int)prompt.size();
            std::cout << "\r\033[" << target << "C" << std::flush;
        }
    }
}

static bool parse_custom_function(const std::string& input, Evaluator& evaluator) {
    size_t ce_pos = input.find(":=");
    if (ce_pos == std::string::npos) return false;
    std::string lhs = trim(input.substr(0, ce_pos));
    std::string rhs = trim(input.substr(ce_pos + 2));
    if (rhs.empty()) return false;

    size_t lp = lhs.find('(');
    size_t rp = lhs.find(')');
    if (lp == std::string::npos || rp == std::string::npos || rp <= lp + 1) return false;

    std::string func_name = trim(lhs.substr(0, lp));
    if (func_name.empty() || is_known_func(func_name)) return false;

    std::string params_str = lhs.substr(lp + 1, rp - lp - 1);
    std::vector<std::string> params;
    std::istringstream ps(params_str);
    std::string p;
    while (std::getline(ps, p, ',')) {
        p = trim(p);
        if (!p.empty()) params.push_back(p);
    }
    if (params.empty()) return false;

    try {
        Lexer lexer(rhs);
        auto tokens = lexer.tokenize();
        Parser parser(tokens);
        ASTPtr body = parser.parse();
        evaluator.set_user_function(func_name, params, body);
        std::cout << "Defined " << func_name << "(";
        for (size_t i = 0; i < params.size(); i++) {
            if (i > 0) std::cout << ",";
            std::cout << params[i];
        }
        std::cout << ") := " << rhs << std::endl;
        return true;
    } catch (...) {
        return false;
    }
}

int main(int argc, char* argv[]) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    GetConsoleMode(hOut, &mode);
    SetConsoleMode(hOut, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#endif

    if (argc >= 3 && std::string(argv[1]) == "-f") {
        std::ifstream file(argv[2]);
        if (!file.is_open()) {
            std::cerr << "Cannot open file: " << argv[2] << std::endl;
            return 1;
        }
        Evaluator evaluator;
        std::string line;
        while (std::getline(file, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;
            try {
                auto parts = split_top_level(line);
                for (size_t pi = 0; pi < parts.size(); pi++) {
                    const std::string& part = parts[pi];
                    if (part.empty()) continue;
                    GraphCmd gcmd = parse_graph_cmd(part);
                    if (gcmd.type != GraphCmd::NONE) continue;
                    if (parse_custom_function(part, evaluator)) continue;
                    size_t eq_pos = part.find('=');
                    if (eq_pos != std::string::npos && eq_pos > 0) {
                        std::string lhs = trim(part.substr(0, eq_pos));
                        std::string rhs_str = trim(part.substr(eq_pos + 1));
                        if (lhs == "x" || lhs == "y" || lhs == "A" || lhs == "B" || lhs == "C" || lhs == "D") {
                            evaluator.set_variable(lhs, evaluator.evaluate(rhs_str));
                            continue;
                        }
                        if (lhs == "VerA" || lhs == "VerB" || lhs == "VerC" || lhs == "VerD") {
                            evaluator.set_vec_variable(lhs, evaluator.evaluate(rhs_str));
                            continue;
                        }
                    }
                    Value result = evaluator.evaluate(part);
                    if (evaluator.latex_mode()) {
                        std::cout << Evaluator::format_latex(result) << std::endl;
                    } else {
                        std::cout << evaluator.format_result(result) << std::endl;
                    }
                }
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
        }
        return 0;
    }

    std::cout << "CCalc v2.0 - Command-line Advanced Calculator" << std::endl;
    std::cout << "Type 'help' for help, 'quit' to exit." << std::endl;
    std::cout << "Up/Down: history, Tab: complete, latex: toggle LaTeX mode" << std::endl;
    std::cout << std::endl;

    Evaluator evaluator;
    std::vector<std::string> history;
    int history_idx = 0;

    while (true) {
        std::string line = readline_custom("CCalc> ", history, history_idx);
        line = trim(line);
        if (line.empty()) continue;

        std::string lower = line;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

        if (lower == "quit" || lower == "exit") {
            std::cout << "Goodbye!" << std::endl;
            break;
        }
        if (lower == "help") {
            print_help();
            continue;
        }
        if (lower == "deg") {
            evaluator.set_angle_mode(Evaluator::DEG);
            std::cout << "Angle mode: degrees" << std::endl;
            continue;
        }
        if (lower == "rad") {
            evaluator.set_angle_mode(Evaluator::RAD);
            std::cout << "Angle mode: radians" << std::endl;
            continue;
        }
        if (lower == "latex") {
            evaluator.set_latex_mode(!evaluator.latex_mode());
            std::cout << "LaTeX mode: " << (evaluator.latex_mode() ? "ON" : "OFF") << std::endl;
            continue;
        }
        if (lower.substr(0, 5) == "base ") {
            try {
                int b = std::stoi(line.substr(5));
                if (b < 2 || b > 36) {
                    std::cout << "Base must be between 2 and 36" << std::endl;
                } else {
                    evaluator.set_output_base(b);
                    std::cout << "Output base set to " << b << std::endl;
                }
            } catch (...) {
                std::cout << "Invalid base value" << std::endl;
            }
            continue;
        }
        if (lower.substr(0, 5) == "prec ") {
            try {
                int p = std::stoi(line.substr(5));
                if (p < 1 || p > 10000) {
                    std::cout << "Precision must be between 1 and 10000" << std::endl;
                } else {
                    g_precision = p;
                    std::cout << "Precision set to " << p << " digits" << std::endl;
                }
            } catch (...) {
                std::cout << "Invalid precision value" << std::endl;
            }
            continue;
        }

        try {
            auto parts = split_top_level(line);
            for (size_t pi = 0; pi < parts.size(); pi++) {
                const std::string& part = parts[pi];
                if (part.empty()) continue;

                if (parse_custom_function(part, evaluator)) continue;

                GraphCmd gcmd = parse_graph_cmd(part);
                if (gcmd.type != GraphCmd::NONE) {
                    std::string cmd;
                    if (gcmd.type == GraphCmd::EXPLICIT) {
                        cmd = "start ccalc_graph.exe \"" + gcmd.expr + "\"";
                        std::cout << "Graph: y = " << gcmd.expr << std::endl;
                    } else if (gcmd.type == GraphCmd::IMPLICIT) {
                        cmd = "start ccalc_graph.exe -i \"" + gcmd.expr + "\"";
                        std::cout << "Graph: " << gcmd.expr << " = 0" << std::endl;
                    } else if (gcmd.type == GraphCmd::PARAMETRIC) {
                        cmd = "start ccalc_graph.exe -p \"" + gcmd.expr + "\" \"" + gcmd.expr2 + "\"";
                        std::cout << "Graph: param(" << gcmd.expr << ", " << gcmd.expr2 << ")" << std::endl;
                    } else if (gcmd.type == GraphCmd::POLAR) {
                        cmd = "start ccalc_graph.exe -l \"" + gcmd.expr + "\"";
                        std::cout << "Graph: r = " << gcmd.expr << std::endl;
                    }
                    std::system(cmd.c_str());
                    continue;
                }

                size_t eq_pos = part.find('=');
                if (eq_pos != std::string::npos && eq_pos > 0) {
                    std::string var_name = trim(part.substr(0, eq_pos));
                    if (is_assignable_var(var_name)) {
                        std::string expr_str = trim(part.substr(eq_pos + 1));
                        if (!expr_str.empty()) {
                            Value result = evaluator.evaluate(expr_str);
                            if (result.is_error()) {
                                std::cout << result.to_string() << std::endl;
                            } else {
                                if (is_vec_var(var_name)) {
                                    evaluator.set_vec_variable(var_name, result);
                                } else {
                                    evaluator.set_variable(var_name, result);
                                }
                                std::cout << var_name << " = "
                                          << Evaluator::format_result(result, evaluator.output_base())
                                          << std::endl;
                            }
                            continue;
                        }
                    }
                }
                Value result = evaluator.evaluate(part);
                if (evaluator.latex_mode()) {
                    std::cout << Evaluator::format_latex(result) << std::endl;
                } else if (result.is_matrix()) {
                    std::cout << Evaluator::format_pretty_matrix(result);
                } else {
                    std::cout << Evaluator::format_result(result, evaluator.output_base())
                              << std::endl;
                }
            }
        } catch (const CalcError& e) {
            std::cout << "Error: " << e.what() << std::endl;
        } catch (const std::exception& e) {
            std::cout << "Error: " << e.what() << std::endl;
        } catch (...) {
            std::cout << "Unknown error" << std::endl;
        }
    }

    return 0;
}
