#include "ccalc.h"
#include <fstream>
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
    std::cout << "  sum(expr, var, from, to)        - Summation" << std::endl;
    std::cout << "  prod(expr, var, from, to)       - Product" << std::endl;
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
    std::cout << "  ans                             - Last answer" << std::endl;
    std::cout << "  help                            - Show this help" << std::endl;
    std::cout << "  quit / exit                     - Exit program" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  sin(pi/2)        => 1" << std::endl;
    std::cout << "  cos(pi/3)        => 1/2" << std::endl;
    std::cout << "  log(2, 8)        => 3" << std::endl;
    std::cout << "  10!              => 3628800" << std::endl;
    std::cout << "  C(10,3)          => 120" << std::endl;
    std::cout << "  sqrt(2)          => sqrt(2) ~= 1.414..." << std::endl;
    std::cout << "  sum(k^2, k, 1, 10)  => 385" << std::endl;
    std::cout << "  x=pi, sin(x/2)   => 1" << std::endl;
    std::cout << "  0xFF             => 255" << std::endl;
    std::cout << "  base 16, 255     => FF" << std::endl;
    std::cout << "  VerA=(1,2,3)     => (1, 2, 3)" << std::endl;
    std::cout << "  vecmod(VerA)      => sqrt(14)" << std::endl;
    std::cout << "  dot(VerA,VerB)    => dot product" << std::endl;
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

int main(int argc, char* argv[]) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
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
                    size_t eq_pos = part.find('=');
                    if (eq_pos != std::string::npos && eq_pos > 0) {
                        std::string lhs = trim(part.substr(0, eq_pos));
                        std::string rhs = trim(part.substr(eq_pos + 1));
                        if (lhs == "x" || lhs == "y" || lhs == "A" || lhs == "B" || lhs == "C" || lhs == "D") {
                            evaluator.set_variable(lhs, evaluator.evaluate(rhs));
                            continue;
                        }
                        if (lhs == "VerA" || lhs == "VerB" || lhs == "VerC" || lhs == "VerD") {
                            evaluator.set_vec_variable(lhs, evaluator.evaluate(rhs));
                            continue;
                        }
                    }
                    Value result = evaluator.evaluate(part);
                    std::cout << evaluator.format_result(result) << std::endl;
                }
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
        }
        return 0;
    }

    std::cout << "CCalc v1.0 - Command-line Advanced Calculator" << std::endl;
    std::cout << "Type 'help' for help, 'quit' to exit." << std::endl;
    std::cout << std::endl;

    Evaluator evaluator;
    std::string line;

    while (true) {
        std::cout << "CCalc> " << std::flush;
        if (!std::getline(std::cin, line)) break;
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
                std::cout << Evaluator::format_result(result, evaluator.output_base())
                          << std::endl;
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
