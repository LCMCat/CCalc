#include "ccalc.h"
#include <complex>
#include <algorithm>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static std::vector<std::complex<double>> solve_quad_d(double a, double b, double c) {
    double D = b * b - 4 * a * c;
    if (D >= 0) {
        double sq = sqrt(D);
        return {{(-b + sq) / (2 * a)}, {(-b - sq) / (2 * a)}};
    }
    double sq = sqrt(-D);
    return {{-b / (2 * a), sq / (2 * a)}, {-b / (2 * a), -sq / (2 * a)}};
}

static double val_to_double(const Value& v) {
    if (v.is_error()) return NAN;
    try { return std::stod(v.to_string()); }
    catch (...) { return NAN; }
}

static std::vector<std::complex<double>> solve_cubic_d(double a, double b, double c, double d);

Value Evaluator::eval_matrix(const std::vector<ASTPtr>& args) {
    if (args.empty()) return Value::make_error("matrix requires arguments");
    Value v = eval_node(args[0]);
    if (v.is_matrix()) return v;
    return Value::make_error("matrix: expected matrix literal or argument");
}

Value Evaluator::eval_det(const Value& m) {
    if (!m.is_matrix()) return Value::make_error("det requires a matrix");
    int n = m.mat_rows();
    if (n != m.mat_cols()) return Value::make_error("det requires a square matrix");
    if (n == 1) return m.mat[0][0];
    if (n == 2) return m.mat[0][0] * m.mat[1][1] - m.mat[0][1] * m.mat[1][0];
    if (n == 3) {
        return m.mat[0][0] * (m.mat[1][1] * m.mat[2][2] - m.mat[1][2] * m.mat[2][1])
             - m.mat[0][1] * (m.mat[1][0] * m.mat[2][2] - m.mat[1][2] * m.mat[2][0])
             + m.mat[0][2] * (m.mat[1][0] * m.mat[2][1] - m.mat[1][1] * m.mat[2][0]);
    }
    Value result(BigRat(0));
    for (int j = 0; j < n; j++) {
        std::vector<std::vector<Value>> sub;
        for (int r = 1; r < n; r++) {
            std::vector<Value> row;
            for (int c = 0; c < n; c++) {
                if (c != j) row.push_back(m.mat[r][c]);
            }
            sub.push_back(row);
        }
        Value minor = eval_det(Value::make_matrix(sub));
        if (minor.is_error()) return minor;
        Value sign(BigRat(j % 2 == 0 ? 1 : -1));
        result = result + sign * m.mat[0][j] * minor;
    }
    return result;
}

Value Evaluator::eval_inv(const Value& m) {
    if (!m.is_matrix()) return Value::make_error("inv requires a matrix");
    int n = m.mat_rows();
    if (n != m.mat_cols()) return Value::make_error("inv requires a square matrix");
    Value d = eval_det(m);
    if (d.is_zero()) return Value::make_error("Matrix is singular (det=0)");
    if (n == 1) {
        Value one(BigRat(1));
        std::vector<std::vector<Value>> r = {{one / m.mat[0][0]}};
        return Value::make_matrix(r);
    }
    std::vector<std::vector<Value>> cofactors;
    for (int i = 0; i < n; i++) {
        std::vector<Value> row;
        for (int j = 0; j < n; j++) {
            std::vector<std::vector<Value>> sub;
            for (int r = 0; r < n; r++) {
                if (r == i) continue;
                std::vector<Value> sr;
                for (int c = 0; c < n; c++) {
                    if (c != j) sr.push_back(m.mat[r][c]);
                }
                sub.push_back(sr);
            }
            Value minor = eval_det(Value::make_matrix(sub));
            Value sign(BigRat((i + j) % 2 == 0 ? 1 : -1));
            row.push_back(sign * minor);
        }
        cofactors.push_back(row);
    }
    std::vector<std::vector<Value>> adj;
    for (int i = 0; i < n; i++) {
        std::vector<Value> row;
        for (int j = 0; j < n; j++) {
            row.push_back(cofactors[j][i] / d);
        }
        adj.push_back(row);
    }
    return Value::make_matrix(adj);
}

Value Evaluator::eval_eigen(const Value& m) {
    if (!m.is_matrix()) return Value::make_error("eigen requires a matrix");
    int n = m.mat_rows();
    if (n != m.mat_cols()) return Value::make_error("eigen requires a square matrix");
    if (n > 4) return Value::make_error("eigen supports up to 4x4 matrices");
    if (n == 1) return Value::make_string("Eigenvalue: " + m.mat[0][0].to_string());
    std::vector<double> coeffs(n + 1, 0);
    if (n == 2) {
        double a = val_to_double(m.mat[0][0]), b = val_to_double(m.mat[0][1]);
        double c = val_to_double(m.mat[1][0]), d = val_to_double(m.mat[1][1]);
        coeffs = {1, -(a + d), a * d - b * c};
    } else if (n == 3) {
        double a = val_to_double(m.mat[0][0]), b = val_to_double(m.mat[0][1]), cc = val_to_double(m.mat[0][2]);
        double d = val_to_double(m.mat[1][0]), e = val_to_double(m.mat[1][1]), f = val_to_double(m.mat[1][2]);
        double g = val_to_double(m.mat[2][0]), h = val_to_double(m.mat[2][1]), k = val_to_double(m.mat[2][2]);
        double tr = a + e + k;
        double det_val = a*(e*k-f*h) - b*(d*k-f*g) + cc*(d*h-e*g);
        double m2 = a*a+b*d+cc*g + d*b+e*e+f*h + g*cc+h*f+k*k;
        coeffs = {1, -tr, (tr*tr - m2)/2, -det_val};
    } else {
        return Value::make_error("eigen for 4x4 not fully supported");
    }
    std::vector<std::complex<double>> roots;
    if (n == 2) roots = solve_quad_d(coeffs[0], coeffs[1], coeffs[2]);
    else if (n == 3) roots = solve_cubic_d(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
    std::string result = "Eigenvalues: ";
    bool first = true;
    for (auto& r : roots) {
        if (!first) result += ", ";
        first = false;
        if (fabs(r.imag()) < 1e-10) {
            char buf[64]; snprintf(buf, sizeof(buf), "%.6g", r.real());
            result += buf;
        } else {
            char buf[128]; snprintf(buf, sizeof(buf), "%.6g%+.6gi", r.real(), r.imag());
            result += buf;
        }
    }
    return Value::make_string(result);
}

Value Evaluator::eval_trace(const Value& m) {
    if (!m.is_matrix()) return Value::make_error("trace requires a matrix");
    if (m.mat_rows() != m.mat_cols()) return Value::make_error("trace requires a square matrix");
    Value sum(BigRat(0));
    for (int i = 0; i < m.mat_rows(); i++)
        sum = sum + m.mat[i][i];
    return sum;
}

Value Evaluator::eval_transpose(const Value& m) {
    if (!m.is_matrix()) return Value::make_error("transpose requires a matrix");
    int r = m.mat_rows(), c = m.mat_cols();
    std::vector<std::vector<Value>> result;
    for (int j = 0; j < c; j++) {
        std::vector<Value> row;
        for (int i = 0; i < r; i++)
            row.push_back(m.mat[i][j]);
        result.push_back(row);
    }
    return Value::make_matrix(result);
}

Value Evaluator::eval_identity(const Value& n) {
    int size = (int)val_to_double(n);
    if (size <= 0 || size > 100) return Value::make_error("identity: invalid size");
    std::vector<std::vector<Value>> result;
    for (int i = 0; i < size; i++) {
        std::vector<Value> row;
        for (int j = 0; j < size; j++)
            row.push_back(Value(BigRat(i == j ? 1 : 0)));
        result.push_back(row);
    }
    return Value::make_matrix(result);
}

Value Evaluator::eval_mean(const std::vector<Value>& args) {
    if (args.empty()) return Value::make_error("mean requires at least 1 argument");
    Value sum(BigRat(0));
    for (auto& v : args) {
        if (v.is_vector()) {
            for (auto& c : v.vec) sum = sum + c;
            sum = sum / Value(BigRat((int)v.vec.size()));
        } else {
            sum = sum + v;
        }
    }
    int total = 0;
    for (auto& v : args) {
        total += v.is_vector() ? (int)v.vec.size() : 1;
    }
    if (total <= 0) return Value::make_error("mean: no values");
    return sum / Value(BigRat(total));
}

Value Evaluator::eval_variance(const std::vector<Value>& args) {
    if (args.size() < 2) return Value::make_error("variance requires at least 2 values");
    std::vector<double> vals;
    for (auto& v : args) {
        if (v.is_vector()) {
            for (auto& c : v.vec) vals.push_back(val_to_double(c));
        } else {
            vals.push_back(val_to_double(v));
        }
    }
    if (vals.size() < 2) return Value::make_error("variance: need at least 2 values");
    double mean = 0;
    for (auto v : vals) mean += v;
    mean /= vals.size();
    double var = 0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    var /= (vals.size() - 1);
    std::ostringstream oss;
    oss << std::setprecision(10) << var;
    return Value(BigRat(BigInt(oss.str()), BigInt(1)));
}

Value Evaluator::eval_stddev(const std::vector<Value>& args) {
    Value v = eval_variance(args);
    if (v.is_error()) return v;
    return eval_sqrt(v);
}

Value Evaluator::eval_median(const std::vector<Value>& args) {
    std::vector<double> vals;
    for (auto& v : args) {
        if (v.is_vector()) {
            for (auto& c : v.vec) vals.push_back(val_to_double(c));
        } else {
            vals.push_back(val_to_double(v));
        }
    }
    if (vals.empty()) return Value::make_error("median: no values");
    std::sort(vals.begin(), vals.end());
    double med;
    size_t n = vals.size();
    if (n % 2 == 1) med = vals[n / 2];
    else med = (vals[n / 2 - 1] + vals[n / 2]) / 2.0;
    std::ostringstream oss;
    oss << std::setprecision(10) << med;
    try { return Value(BigRat(oss.str())); }
    catch (...) { return Value(BigFloat(med)); }
}

static std::string ast_to_string(ASTPtr node) {
    if (!node) return "";
    switch (node->type) {
    case ASTNode::NUMBER: return node->number.to_string();
    case ASTNode::CONSTANT: return node->name;
    case ASTNode::VARIABLE: return node->name;
    case ASTNode::BINOP: {
        std::string l = ast_to_string(node->left);
        std::string r = ast_to_string(node->right);
        if (node->op == '+') return l + "+" + r;
        if (node->op == '-') return l + "-" + r;
        if (node->op == '*') return "(" + l + ")*(" + r + ")";
        if (node->op == '/') return "(" + l + ")/(" + r + ")";
        if (node->op == '^') return "(" + l + ")^(" + r + ")";
        return l + "?" + r;
    }
    case ASTNode::UNARYOP:
        return "-(" + ast_to_string(node->left) + ")";
    case ASTNode::FUNCTION: {
        std::string r = node->name + "(";
        for (size_t i = 0; i < node->args.size(); i++) {
            if (i > 0) r += ",";
            r += ast_to_string(node->args[i]);
        }
        return r + ")";
    }
    default: return "?";
    }
}

static bool is_poly_in_var(ASTPtr node, const std::string& var) {
    if (!node) return true;
    switch (node->type) {
    case ASTNode::NUMBER: case ASTNode::CONSTANT: return true;
    case ASTNode::VARIABLE: return node->name != var || true;
    case ASTNode::BINOP: return is_poly_in_var(node->left, var) && is_poly_in_var(node->right, var);
    case ASTNode::UNARYOP: return is_poly_in_var(node->left, var);
    case ASTNode::FUNCTION: {
        for (auto& a : node->args) if (!is_poly_in_var(a, var)) return false;
        return true;
    }
    default: return false;
    }
}

Value Evaluator::eval_simplify(ASTPtr node) {
    Value result = eval_node(node);
    if (result.is_error()) return result;
    if (!node) return result;

    if (node->type == ASTNode::BINOP && node->op == '/') {
        Value num_val = eval_node(node->left);
        Value den_val = eval_node(node->right);
        if (!num_val.is_error() && !den_val.is_error()) {
            Value simplified = num_val / den_val;
            if (!simplified.is_error()) {
                return Value::make_string(ast_to_string(node->left) + "/" + ast_to_string(node->right)
                    + " = " + simplified.to_string());
            }
        }
    }

    return Value::make_string(ast_to_string(node) + " = " + result.to_string());
}

static std::vector<std::complex<double>> solve_cubic_d(double a, double b,
                                                        double c, double d) {
    if (fabs(a) < 1e-15) return solve_quad_d(b, c, d);
    double p = (3 * a * c - b * b) / (3 * a * a);
    double q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);
    double sh = -b / (3 * a);
    double Q = p * p * p / 27 + q * q / 4;
    if (Q > 1e-14) {
        double sq = sqrt(Q);
        double al = cbrt(-q / 2 + sq);
        double be = cbrt(-q / 2 - sq);
        std::complex<double> w(-0.5, sqrt(3.0) / 2);
        return {al + be + sh, w * al + std::conj(w) * be + sh,
                std::conj(w) * al + w * be + sh};
    }
    if (Q < -1e-14) {
        double mp = -p;
        double m = 2 * sqrt(mp / 3);
        double arg = -3 * q * sqrt(3.0) / (2 * pow(mp, 1.5));
        if (arg > 1.0) arg = 1.0;
        if (arg < -1.0) arg = -1.0;
        double th = acos(arg) / 3;
        return {m * cos(th) + sh, m * cos(th - 2 * M_PI / 3) + sh,
                m * cos(th - 4 * M_PI / 3) + sh};
    }
    double u = cbrt(-q / 2);
    return {2 * u + sh, -u + sh, -u + sh};
}

static std::vector<std::complex<double>> solve_quartic_d(double a, double b,
                                                          double c, double d,
                                                          double e) {
    if (fabs(a) < 1e-15) return solve_cubic_d(b, c, d, e);
    double sh = -b / (4 * a);
    double aa = a;
    b /= aa; c /= aa; d /= aa; e /= aa;
    double p = c - 3 * b * b / 8;
    double q = b * b * b / 8 - b * c / 2 + d;
    double r = -3 * b * b * b * b / 256 + b * b * c / 16 - b * d / 4 + e;
    if (fabs(q) < 1e-14) {
        auto zr = solve_quad_d(1, p, r);
        std::vector<std::complex<double>> roots;
        for (auto& z : zr) {
            roots.push_back(sqrt(z));
            roots.push_back(-sqrt(z));
        }
        for (auto& rt : roots) rt += sh;
        return roots;
    }
    auto sr = solve_cubic_d(8, -4 * p, -8 * r, 4 * p * r - q * q);
    double s = 0;
    double mi = 1e30;
    for (auto& rt : sr) {
        if (fabs(rt.imag()) < mi) {
            mi = fabs(rt.imag());
            s = rt.real();
        }
    }
    double t2 = 2 * s - p;
    std::vector<std::complex<double>> roots;
    if (t2 >= 0) {
        double t = sqrt(t2);
        double u = -q / (2 * t);
        auto q1 = solve_quad_d(1, t, s + u);
        auto q2 = solve_quad_d(1, -t, s - u);
        for (auto& r : q1) roots.push_back(r + sh);
        for (auto& r : q2) roots.push_back(r + sh);
    } else {
        std::complex<double> tc = sqrt(std::complex<double>(t2));
        std::complex<double> uc = -q / (2.0 * tc);
        std::complex<double> sqD1 = sqrt(tc * tc - 4.0 * (s + uc));
        std::complex<double> sqD2 = sqrt(tc * tc - 4.0 * (s - uc));
        roots.push_back((-tc + sqD1) / 2.0 + sh);
        roots.push_back((-tc - sqD1) / 2.0 + sh);
        roots.push_back((tc + sqD2) / 2.0 + sh);
        roots.push_back((tc - sqD2) / 2.0 + sh);
    }
    return roots;
}

Value Evaluator::evaluate(ASTPtr node) {
    return eval_node(node);
}

Value Evaluator::evaluate(const std::string& input) {
    Lexer lexer(input);
    auto tokens = lexer.tokenize();
    Parser parser(tokens);
    auto ast = parser.parse();
    Value result = evaluate(ast);
    last_ans_ = result;
    return result;
}

void Evaluator::set_variable(const std::string& name, const Value& v) {
    variables_[name] = v;
}

Value Evaluator::get_variable(const std::string& name) const {
    if (name == "ans") return last_ans_;
    auto it = variables_.find(name);
    if (it != variables_.end()) return it->second;
    auto vit = vec_variables_.find(name);
    if (vit != vec_variables_.end()) return vit->second;
    throw CalcError("Undefined variable: " + name);
}

void Evaluator::set_vec_variable(const std::string& name, const Value& v) {
    vec_variables_[name] = v;
}

Value Evaluator::get_vec_variable(const std::string& name) const {
    auto it = vec_variables_.find(name);
    if (it != vec_variables_.end()) return it->second;
    return Value::make_vector({});
}

bool Evaluator::has_vec_variable(const std::string& name) const {
    return vec_variables_.find(name) != vec_variables_.end();
}

Value Evaluator::eval_node(ASTPtr node) {
    switch (node->type) {
    case ASTNode::NUMBER:
        return Value(node->number);
    case ASTNode::CONSTANT:
        if (node->name == "pi") {
            return Value(SurdsExpr(BigRat(1), BigInt(-1)));
        }
        if (node->name == "e") {
            return Value(SurdsExpr(BigRat(1), BigInt(-2)));
        }
        throw CalcError("Unknown constant: " + node->name);
    case ASTNode::VARIABLE:
        return get_variable(node->name);
    case ASTNode::BINOP:
        return eval_binop(node);
    case ASTNode::UNARYOP:
        return eval_unaryop(node);
    case ASTNode::FUNCTION:
        return eval_function(node);
    case ASTNode::FACTORIAL:
        return eval_factorial(node);
    case ASTNode::VEC_LITERAL: {
        std::vector<Value> components;
        for (auto& arg : node->args) {
            Value v = eval_node(arg);
            if (v.is_error()) return v;
            components.push_back(v);
        }
        return Value::make_vector(components);
    }
    case ASTNode::MAT_LITERAL: {
        std::vector<std::vector<Value>> rows;
        for (auto& row : node->mat_rows) {
            std::vector<Value> r;
            for (auto& arg : row) {
                Value v = eval_node(arg);
                if (v.is_error()) return v;
                r.push_back(v);
            }
            rows.push_back(r);
        }
        return Value::make_matrix(rows);
    }
    }
    return Value::make_error("Unknown node type");
}

Value Evaluator::eval_binop(ASTPtr node) {
    Value left = eval_node(node->left);
    Value right = eval_node(node->right);
    if (left.is_error()) return left;
    if (right.is_error()) return right;
    switch (node->op) {
    case '+': return left + right;
    case '-': return left - right;
    case '*': return left * right;
    case '/': return left / right;
    case '^': {
        if (left.is_zero() && right.is_zero())
            return Value::make_error("0^0 is undefined");
        if (right.is_zero()) return Value(BigRat(1));
        if (left.is_zero()) return Value(BigRat(0));
        if (left.type == Value::SURDS && right.type == Value::SURDS) {
            if (right.surds.is_rational()) {
                BigRat exp_r = right.surds.to_rational();
                if (exp_r.denominator() == BigInt(1)) {
                    int64_t n = exp_r.numerator().to_int64();
                    if (left.surds.is_rational()) {
                        BigRat base = left.surds.to_rational();
                        return Value(BigRat::pow(base, n));
                    }
                    bool left_is_e = false;
                    if (left.surds.terms.size() == 1 &&
                        left.surds.terms[0].radicand == BigInt(-2)) {
                        left_is_e = true;
                    }
                    if (n == 1) return Value(left.surds);
                    if (left_is_e) {
                        return Value(BigFloat::exp(BigFloat(n)));
                    }
                    if (n > 1) {
                        SurdsExpr result(BigRat(1));
                        for (int64_t i = 0; i < n; i++) {
                            SurdsExpr prev = result;
                            result = result * left.surds;
                            if (result.is_zero() && !prev.is_zero()) {
                                return Value(BigFloat::pow_val(left.to_float(), right.to_float()));
                            }
                        }
                        if (result.is_rational()) return Value(result.to_rational());
                        return Value(result);
                    }
                }
            }
        }
        return Value(BigFloat::pow_val(left.to_float(), right.to_float()));
    }
    case '%': {
        if (left.is_rational() && right.is_rational()) {
            BigRat lr = left.to_rational(), rr = right.to_rational();
            if (rr.denominator() == BigInt(1) && lr.denominator() == BigInt(1)) {
                BigInt ln = lr.numerator(), rn = rr.numerator();
                if (rn.is_zero()) return Value::make_error("Modulo by zero");
                return Value(ln % rn);
            }
        }
        BigFloat lf = left.to_float(), rf = right.to_float();
        BigFloat q = lf / rf;
        std::string qs = q.to_string(30);
        size_t dot = qs.find('.');
        if (dot != std::string::npos) qs = qs.substr(0, dot);
        BigFloat qi(qs);
        return Value(lf - qi * rf);
    }
    case '<': case '>': case '=': case 'L': case 'G': {
        BigFloat lf = left.to_float();
        BigFloat rf = right.to_float();
        bool result = false;
        switch (node->op) {
        case '<': result = lf < rf; break;
        case '>': result = lf > rf; break;
        case '=': result = lf == rf; break;
        case 'L': result = lf <= rf; break;
        case 'G': result = lf >= rf; break;
        }
        return Value(BigRat(result ? 1 : 0));
    }
    }
    return Value::make_error("Unknown operator");
}

Value Evaluator::eval_unaryop(ASTPtr node) {
    Value operand = eval_node(node->left);
    if (operand.is_error()) return operand;
    if (node->op == '-') return -operand;
    return operand;
}

Value Evaluator::eval_factorial(ASTPtr node) {
    return eval_fact(eval_node(node->left));
}

Value Evaluator::to_radians(const Value& v) const {
    if (angle_mode_ == DEG) {
        if (v.is_rational()) {
            BigRat r = v.to_rational();
            BigRat pi_coeff = r / BigRat(180);
            return Value(SurdsExpr(pi_coeff, BigInt(-1)));
        }
        return Value(v.to_float() * BigFloat::pi() / BigFloat(180));
    }
    return v;
}

Value Evaluator::from_radians(const Value& v) const {
    if (angle_mode_ == DEG) {
        if (v.type == Value::SURDS && v.surds.terms.size() == 1 &&
            v.surds.terms[0].radicand == BigInt(-1)) {
            BigRat pi_coeff = v.surds.terms[0].coeff;
            BigRat deg = pi_coeff * BigRat(180);
            if (deg.denominator() == BigInt(1)) return Value(deg);
            return Value(deg);
        }
        return Value(v.to_float() * BigFloat(180) / BigFloat::pi());
    }
    return v;
}

Value Evaluator::try_exact_trig(const BigRat& pi_coeff, int func) {
    BigRat q = pi_coeff;
    while (q < BigRat(0)) q = q + BigRat(2);
    while (q >= BigRat(2)) q = q - BigRat(2);
    if (q == BigRat(0)) {
        if (func == 0) return Value(BigRat(0));
        if (func == 1) return Value(BigRat(1));
        if (func == 2) return Value(BigRat(0));
    }
    if (q == BigRat(1, 6)) {
        if (func == 0) return Value(BigRat(1, 2));
        if (func == 1) return Value(SurdsExpr(BigRat(1, 2), BigInt(3)));
        if (func == 2) return Value(SurdsExpr(BigRat(1, 3), BigInt(3)));
    }
    if (q == BigRat(1, 4)) {
        if (func == 0) return Value(SurdsExpr(BigRat(1, 2), BigInt(2)));
        if (func == 1) return Value(SurdsExpr(BigRat(1, 2), BigInt(2)));
        if (func == 2) return Value(BigRat(1));
    }
    if (q == BigRat(1, 3)) {
        if (func == 0) return Value(SurdsExpr(BigRat(1, 2), BigInt(3)));
        if (func == 1) return Value(BigRat(1, 2));
        if (func == 2) return Value(SurdsExpr(BigRat(1), BigInt(3)));
    }
    if (q == BigRat(1, 2)) {
        if (func == 0) return Value(BigRat(1));
        if (func == 1) return Value(BigRat(0));
        if (func == 2) return Value::make_error("tan(pi/2) is undefined");
    }
    if (q == BigRat(2, 3)) {
        if (func == 0) return Value(SurdsExpr(BigRat(1, 2), BigInt(3)));
        if (func == 1) return Value(BigRat(-1, 2));
        if (func == 2) return Value(SurdsExpr(BigRat(-1), BigInt(3)));
    }
    if (q == BigRat(3, 4)) {
        if (func == 0) return Value(SurdsExpr(BigRat(1, 2), BigInt(2)));
        if (func == 1) return Value(SurdsExpr(BigRat(-1, 2), BigInt(2)));
        if (func == 2) return Value(BigRat(-1));
    }
    if (q == BigRat(5, 6)) {
        if (func == 0) return Value(BigRat(1, 2));
        if (func == 1) return Value(SurdsExpr(BigRat(-1, 2), BigInt(3)));
        if (func == 2) return Value(SurdsExpr(BigRat(-1, 3), BigInt(3)));
    }
    if (q == BigRat(1)) {
        if (func == 0) return Value(BigRat(0));
        if (func == 1) return Value(BigRat(-1));
        if (func == 2) return Value(BigRat(0));
    }
    if (q > BigRat(1)) {
        BigRat q2 = BigRat(2) - q;
        Value v = try_exact_trig(q2, func);
        if (v.is_error()) return v;
        if (func == 0) return -v;
        if (func == 1) return v;
        if (func == 2) return -v;
    }
    if (q == BigRat(1, 12)) {
        SurdsExpr sin_r;
        sin_r.terms = {SurdsTerm(BigRat(1, 4), BigInt(6)), SurdsTerm(BigRat(-1, 4), BigInt(2))};
        SurdsExpr cos_r;
        cos_r.terms = {SurdsTerm(BigRat(1, 4), BigInt(6)), SurdsTerm(BigRat(1, 4), BigInt(2))};
        if (func == 0) return Value(sin_r);
        if (func == 1) return Value(cos_r);
        if (func == 2) return Value(sin_r) / Value(cos_r);
    }
    if (q == BigRat(5, 12)) {
        SurdsExpr sin_r;
        sin_r.terms = {SurdsTerm(BigRat(1, 4), BigInt(6)), SurdsTerm(BigRat(1, 4), BigInt(2))};
        SurdsExpr cos_r;
        cos_r.terms = {SurdsTerm(BigRat(1, 4), BigInt(6)), SurdsTerm(BigRat(-1, 4), BigInt(2))};
        if (func == 0) return Value(sin_r);
        if (func == 1) return Value(cos_r);
        if (func == 2) return Value(sin_r) / Value(cos_r);
    }
    if (q == BigRat(7, 12)) {
        SurdsExpr sin_r;
        sin_r.terms = {SurdsTerm(BigRat(1, 4), BigInt(6)), SurdsTerm(BigRat(1, 4), BigInt(2))};
        SurdsExpr cos_r;
        cos_r.terms = {SurdsTerm(BigRat(-1, 4), BigInt(6)), SurdsTerm(BigRat(1, 4), BigInt(2))};
        if (func == 0) return Value(sin_r);
        if (func == 1) return Value(cos_r);
        if (func == 2) return Value(sin_r) / Value(cos_r);
    }
    if (q == BigRat(11, 12)) {
        SurdsExpr sin_r;
        sin_r.terms = {SurdsTerm(BigRat(1, 4), BigInt(6)), SurdsTerm(BigRat(-1, 4), BigInt(2))};
        SurdsExpr cos_r;
        cos_r.terms = {SurdsTerm(BigRat(-1, 4), BigInt(6)), SurdsTerm(BigRat(-1, 4), BigInt(2))};
        if (func == 0) return Value(sin_r);
        if (func == 1) return Value(cos_r);
        if (func == 2) return Value(sin_r) / Value(cos_r);
    }
    return Value();
}

Value Evaluator::eval_sin(const Value& v) {
    Value angle = to_radians(v);
    if (angle.type == Value::SURDS && angle.surds.is_rational()) {
        BigRat r = angle.surds.to_rational();
        if (r.is_zero()) return Value(BigRat(0));
    }
    if (angle.type == Value::SURDS && angle.surds.terms.size() == 1) {
        auto& t = angle.surds.terms[0];
        if (t.radicand == BigInt(-1)) {
            BigRat pi_coeff = t.coeff;
            Value exact = try_exact_trig(pi_coeff, 0);
            if (!exact.is_error()) return exact;
        }
    }
    return Value(BigFloat::sin_val(angle.to_float()));
}

Value Evaluator::eval_cos(const Value& v) {
    Value angle = to_radians(v);
    if (angle.type == Value::SURDS && angle.surds.is_rational()) {
        BigRat r = angle.surds.to_rational();
        if (r.is_zero()) return Value(BigRat(1));
    }
    if (angle.type == Value::SURDS && angle.surds.terms.size() == 1) {
        auto& t = angle.surds.terms[0];
        if (t.radicand == BigInt(-1)) {
            BigRat pi_coeff = t.coeff;
            Value exact = try_exact_trig(pi_coeff, 1);
            if (!exact.is_error()) return exact;
        }
    }
    return Value(BigFloat::cos_val(angle.to_float()));
}

Value Evaluator::eval_tan(const Value& v) {
    Value angle = to_radians(v);
    if (angle.type == Value::SURDS && angle.surds.is_rational()) {
        BigRat r = angle.surds.to_rational();
        if (r.is_zero()) return Value(BigRat(0));
    }
    if (angle.type == Value::SURDS && angle.surds.terms.size() == 1) {
        auto& t = angle.surds.terms[0];
        if (t.radicand == BigInt(-1)) {
            BigRat pi_coeff = t.coeff;
            Value exact = try_exact_trig(pi_coeff, 2);
            if (!exact.is_error()) return exact;
        }
    }
    return Value(BigFloat::tan_val(angle.to_float()));
}

Value Evaluator::eval_asin(const Value& v) {
    if (v.type == Value::SURDS && v.surds.is_rational()) {
        BigRat r = v.surds.to_rational();
        if (r == BigRat(0)) return Value(BigRat(0));
        if (r == BigRat(1)) return Value(SurdsExpr(BigRat(1, 2), BigInt(-1)));
        if (r == BigRat(-1)) return Value(SurdsExpr(BigRat(-1, 2), BigInt(-1)));
        if (r == BigRat(1, 2)) return Value(SurdsExpr(BigRat(1, 6), BigInt(-1)));
        if (r == BigRat(-1, 2)) return Value(SurdsExpr(BigRat(-1, 6), BigInt(-1)));
    }
    Value result(Value(BigFloat::asin_val(v.to_float())));
    return from_radians(result);
}

Value Evaluator::eval_acos(const Value& v) {
    if (v.type == Value::SURDS && v.surds.is_rational()) {
        BigRat r = v.surds.to_rational();
        if (r == BigRat(1)) return Value(BigRat(0));
        if (r == BigRat(-1)) return Value(SurdsExpr(BigRat(1), BigInt(-1)));
        if (r == BigRat(0)) return Value(SurdsExpr(BigRat(1, 2), BigInt(-1)));
        if (r == BigRat(1, 2)) return Value(SurdsExpr(BigRat(1, 3), BigInt(-1)));
        if (r == BigRat(-1, 2)) return Value(SurdsExpr(BigRat(2, 3), BigInt(-1)));
    }
    Value result(Value(BigFloat::acos_val(v.to_float())));
    return from_radians(result);
}

Value Evaluator::eval_atan(const Value& v) {
    if (v.type == Value::SURDS && v.surds.is_rational()) {
        BigRat r = v.surds.to_rational();
        if (r == BigRat(0)) return Value(BigRat(0));
        if (r == BigRat(1)) return Value(SurdsExpr(BigRat(1, 4), BigInt(-1)));
        if (r == BigRat(-1)) return Value(SurdsExpr(BigRat(-1, 4), BigInt(-1)));
    }
    Value result(Value(BigFloat::atan_val(v.to_float())));
    return from_radians(result);
}

Value Evaluator::eval_sinh(const Value& v) {
    return Value(BigFloat::sinh_val(v.to_float()));
}

Value Evaluator::eval_cosh(const Value& v) {
    return Value(BigFloat::cosh_val(v.to_float()));
}

Value Evaluator::eval_tanh(const Value& v) {
    return Value(BigFloat::tanh_val(v.to_float()));
}

Value Evaluator::eval_sqrt(const Value& v) {
    if (v.is_zero()) return Value(BigRat(0));
    if (v.type == Value::SURDS) {
        if (v.surds.is_rational()) {
            BigRat r = v.surds.to_rational();
            if (r < BigRat(0)) {
                Value neg_v = Value(-r);
                Value sq = eval_sqrt(neg_v);
                return Value::make_complex(Value(BigRat(0)), sq);
            }
            if (r.numerator().is_zero()) return Value(BigRat(0));
            BigInt num = r.numerator().abs();
            BigInt den = r.denominator();
            BigInt a_num(1), a_den(1);
            BigInt r_num = num, r_den = den;
            for (int64_t i = 2; i * i <= 1000000; i++) {
                BigInt ii(i);
                BigInt ii2 = ii * ii;
                while (true) {
                    auto dm = r_num.divmod(ii2);
                    if (dm.second.is_zero()) {
                        r_num = dm.first;
                        a_num = a_num * ii;
                    } else break;
                }
                while (true) {
                    auto dm = r_den.divmod(ii2);
                    if (dm.second.is_zero()) {
                        r_den = dm.first;
                        a_den = a_den * ii;
                    } else break;
                }
            }
            BigRat coeff(a_num, a_den);
            if (r_num == BigInt(1) && r_den == BigInt(1)) {
                return Value(coeff);
            }
            if (r_den == BigInt(1)) {
                if (coeff == BigRat(1)) return Value(SurdsExpr(BigRat(1), r_num));
                return Value(SurdsExpr(coeff, r_num));
            }
            BigInt combined = r_num * r_den;
            BigRat new_coeff = coeff * BigRat(BigInt(1), r_den);
            if (combined == BigInt(1)) return Value(new_coeff);
            return Value(SurdsExpr(new_coeff, combined));
        }
    }
    if (v.type == Value::COMPLEX || v.to_float() < BigFloat(0)) {
        Value neg_v = -v;
        Value sq = eval_sqrt(neg_v);
        return Value::make_complex(Value(BigRat(0)), sq);
    }
    return Value(BigFloat::sqrt_val(v.to_float()));
}

Value Evaluator::eval_cbrt(const Value& v) {
    if (v.is_zero()) return Value(BigRat(0));
    return Value(BigFloat::cbrt_val(v.to_float()));
}

Value Evaluator::eval_nrt(const Value& v, const Value& n) {
    if (!n.is_rational())
        return Value::make_error("n-th root requires integer n");
    BigRat nr = n.to_rational();
    if (nr.denominator() != BigInt(1))
        return Value::make_error("n-th root requires integer n");
    int64_t nv = nr.numerator().to_int64();
    return Value(BigFloat::nrt_val(v.to_float(), nv));
}

Value Evaluator::eval_abs(const Value& v) {
    return v.magnitude();
}

Value Evaluator::eval_ln(const Value& v) {
    if (v.is_zero()) return Value::make_error("ln(0) is undefined");
    if (v.type == Value::SURDS) {
        if (v.surds.is_rational()) {
            BigRat r = v.surds.to_rational();
            if (r == BigRat(1)) return Value(BigRat(0));
            if (r == BigRat(-1)) {
                return Value::make_complex(Value(BigRat(0)), Value(BigFloat::pi()));
            }
        }
        if (v.surds.terms.size() == 1) {
            auto& t = v.surds.terms[0];
            if (t.radicand == BigInt(-2)) {
                if (t.coeff == BigRat(1)) return Value(BigRat(1));
                if (t.coeff == BigRat(-1)) return Value(BigRat(-1));
            }
        }
    }
    if (v.type == Value::COMPLEX) {
        Value mag = v.magnitude();
        Value arg = v.argument();
        return Value::make_complex(eval_ln(mag), arg);
    }
    if (v.to_float() < BigFloat(0)) {
        Value pos = -v;
        Value ln_pos = eval_ln(pos);
        return Value::make_complex(ln_pos, Value(BigFloat::pi()));
    }
    return Value(BigFloat::ln(v.to_float()));
}

Value Evaluator::eval_lg(const Value& v) {
    if (v.is_zero()) return Value::make_error("lg(0) is undefined");
    if (v.type == Value::SURDS && v.surds.is_rational()) {
        BigRat r = v.surds.to_rational();
        if (r == BigRat(1)) return Value(BigRat(0));
        if (r == BigRat(10)) return Value(BigRat(1));
        if (r == BigRat(100)) return Value(BigRat(2));
        if (r == BigRat(1000)) return Value(BigRat(3));
        if (r.numerator() == BigInt(10)) {
            BigInt d = r.denominator();
            std::string ds = d.to_string();
            bool is_pow_10 = (ds[0] == '1');
            for (size_t i = 1; i < ds.size(); i++) if (ds[i] != '0') is_pow_10 = false;
            if (is_pow_10) return Value(BigRat(1) - BigRat((int64_t)ds.size() - 1));
        }
    }
    return Value(BigFloat::log10(v.to_float()));
}

Value Evaluator::eval_log(const Value& base, const Value& v) {
    if (v.is_zero()) return Value::make_error("log(0) is undefined");
    if (v.type == Value::SURDS && v.surds.terms.size() == 1 &&
        v.surds.terms[0].radicand == BigInt(-2) &&
        base.type == Value::SURDS && base.surds.terms.size() == 1 &&
        base.surds.terms[0].radicand == BigInt(-2)) {
        return Value(BigRat(1));
    }
    if (base.type == Value::SURDS && base.surds.is_rational() &&
        v.type == Value::SURDS && v.surds.is_rational()) {
        BigRat b = base.surds.to_rational();
        BigRat x = v.surds.to_rational();
        if (x == BigRat(1)) return Value(BigRat(0));
        if (x == b) return Value(BigRat(1));
        if (b == BigRat(10)) return eval_lg(v);
        if (b == BigRat(2)) {
            if (x == BigRat(2)) return Value(BigRat(1));
            if (x == BigRat(4)) return Value(BigRat(2));
            if (x == BigRat(8)) return Value(BigRat(3));
            if (x == BigRat(16)) return Value(BigRat(4));
            if (x == BigRat(1, 2)) return Value(BigRat(-1));
            if (x == BigRat(1, 4)) return Value(BigRat(-2));
        }
    }
    return Value(BigFloat::log_base(v.to_float(), base.to_float()));
}

Value Evaluator::eval_fact(const Value& v) {
    if (!v.is_rational())
        return Value::make_error("Factorial requires non-negative integer");
    BigRat r = v.to_rational();
    if (r.denominator() != BigInt(1) || r < BigRat(0))
        return Value::make_error("Factorial requires non-negative integer");
    int64_t n = r.numerator().to_int64();
    return Value(BigInt::factorial(n));
}

Value Evaluator::eval_perm(const Value& n, const Value& k) {
    if (!n.is_rational() || !k.is_rational())
        return Value::make_error("Permutation requires non-negative integers");
    BigRat nr = n.to_rational(), kr = k.to_rational();
    if (nr.denominator() != BigInt(1) || kr.denominator() != BigInt(1))
        return Value::make_error("Permutation requires non-negative integers");
    int64_t nv = nr.numerator().to_int64();
    int64_t kv = kr.numerator().to_int64();
    if (kv > nv || kv < 0)
        return Value::make_error("Invalid permutation arguments");
    BigInt result(BigInt(1));
    for (int64_t i = nv; i > nv - kv; i--) result = result * BigInt(i);
    return Value(result);
}

Value Evaluator::eval_comb(const Value& n, const Value& k) {
    if (!n.is_rational() || !k.is_rational())
        return Value::make_error("Combination requires non-negative integers");
    BigRat nr = n.to_rational(), kr = k.to_rational();
    if (nr.denominator() != BigInt(1) || kr.denominator() != BigInt(1))
        return Value::make_error("Combination requires non-negative integers");
    int64_t nv = nr.numerator().to_int64();
    int64_t kv = kr.numerator().to_int64();
    if (kv > nv || kv < 0)
        return Value::make_error("Invalid combination arguments");
    if (kv > nv - kv) kv = nv - kv;
    BigInt result(BigInt(1));
    for (int64_t i = 0; i < kv; i++) {
        result = result * BigInt(nv - i) / BigInt(i + 1);
    }
    return Value(result);
}

Value Evaluator::eval_gcd(const Value& a, const Value& b) {
    if (!a.is_rational() || !b.is_rational())
        return Value::make_error("GCD requires integers");
    BigRat ar = a.to_rational(), br = b.to_rational();
    if (ar.denominator() != BigInt(1) || br.denominator() != BigInt(1))
        return Value::make_error("GCD requires integers");
    return Value(BigInt::gcd(ar.numerator(), br.numerator()));
}

Value Evaluator::eval_lcm(const Value& a, const Value& b) {
    if (!a.is_rational() || !b.is_rational())
        return Value::make_error("LCM requires integers");
    BigRat ar = a.to_rational(), br = b.to_rational();
    if (ar.denominator() != BigInt(1) || br.denominator() != BigInt(1))
        return Value::make_error("LCM requires integers");
    return Value(BigInt::lcm(ar.numerator(), br.numerator()));
}

Value Evaluator::eval_factor(const Value& v) {
    if (!v.is_rational())
        return Value::make_error("Factorization requires positive integer");
    BigRat r = v.to_rational();
    if (r.denominator() != BigInt(1) || r < BigRat(0))
        return Value::make_error("Factorization requires positive integer");
    BigInt n = r.numerator();
    if (n == BigInt(1)) return Value(BigRat(1));
    std::map<int64_t, int> factors;
    int64_t nv = n.to_int64();
    if (nv == 0) {
        for (int64_t d = 2; d * d <= 1000000000LL && n > BigInt(1); d++) {
            int count = 0;
            while (n.divmod(BigInt(d)).second.is_zero()) {
                n = n / BigInt(d);
                count++;
            }
            if (count > 0) factors[d] = count;
        }
    } else {
        for (int64_t d = 2; d * d <= nv; d++) {
            int count = 0;
            while (nv % d == 0) { nv /= d; count++; }
            if (count > 0) factors[d] = count;
        }
        if (nv > 1) factors[nv] = 1;
    }
    if (n > BigInt(1) && nv == 0) {
        for (int64_t d = 1000000001LL; BigInt(d) * BigInt(d) <= n; d += 2) {
            int count = 0;
            while (n.divmod(BigInt(d)).second.is_zero()) {
                n = n / BigInt(d);
                count++;
            }
            if (count > 0) factors[d] = count;
        }
        if (n > BigInt(1)) factors[n.to_int64() > 0 ? n.to_int64() : 0] = 1;
    }
    std::string result;
    for (auto& [p, e] : factors) {
        if (!result.empty()) result += " * ";
        result += std::to_string(p);
        if (e > 1) result += "^" + std::to_string(e);
    }
    return Value::make_string(result);
}

Value Evaluator::substitute(ASTPtr node, const std::string& var, const Value& val) {
    if (!node) return Value::make_error("Null node");
    switch (node->type) {
    case ASTNode::NUMBER:
        return Value(node->number);
    case ASTNode::CONSTANT:
        if (node->name == "pi") return Value(SurdsExpr(BigRat(1), BigInt(-1)));
        if (node->name == "e") return Value(SurdsExpr(BigRat(1), BigInt(-2)));
        throw CalcError("Unknown constant: " + node->name);
    case ASTNode::VARIABLE:
        if (node->name == var) return val;
        return get_variable(node->name);
    case ASTNode::BINOP: {
        Value l = substitute(node->left, var, val);
        Value r = substitute(node->right, var, val);
        if (l.is_error()) return l;
        if (r.is_error()) return r;
        switch (node->op) {
        case '+': return l + r;
        case '-': return l - r;
        case '*': return l * r;
        case '/': return l / r;
        case '^': {
            if (l.is_zero() && r.is_zero()) return Value::make_error("0^0");
            if (r.is_zero()) return Value(BigRat(1));
            if (l.is_zero()) return Value(BigRat(0));
            if (r.is_rational()) {
                BigRat exp_r = r.to_rational();
                if (exp_r.denominator() == BigInt(1)) {
                    int64_t n = exp_r.numerator().to_int64();
                    if (l.is_rational()) {
                        return Value(BigRat::pow(l.to_rational(), n));
                    }
                    if (l.type == Value::FLOAT) {
                        if (n >= 0 && n <= 1000) {
                            BigFloat result(1);
                            BigFloat base = l.to_float();
                            for (int64_t i = 0; i < n; i++) result = result * base;
                            return Value(result);
                        }
                    }
                    if (l.type == Value::SURDS) {
                        bool left_is_e = (l.surds.terms.size() == 1 &&
                            l.surds.terms[0].radicand == BigInt(-2));
                        if (n == 1) return Value(l.surds);
                        if (left_is_e) return Value(BigFloat::exp(BigFloat(n)));
                        if (n > 1) {
                            SurdsExpr result(BigRat(1));
                            for (int64_t i = 0; i < n; i++) {
                                SurdsExpr prev = result;
                                result = result * l.surds;
                                if (result.is_zero() && !prev.is_zero())
                                    return Value(BigFloat::pow_val(l.to_float(), r.to_float()));
                            }
                            if (result.is_rational()) return Value(result.to_rational());
                            return Value(result);
                        }
                    }
                }
            }
            return Value(BigFloat::pow_val(l.to_float(), r.to_float()));
        }
        case '%': {
            if (r.is_zero()) return Value::make_error("Division by zero");
            if (l.type == Value::SURDS && r.type == Value::SURDS &&
                l.surds.is_rational() && r.surds.is_rational()) {
                BigRat lr = l.surds.to_rational();
                BigRat rr = r.surds.to_rational();
                if (lr.denominator() == BigInt(1) && rr.denominator() == BigInt(1)) {
                    auto dm = lr.numerator().divmod(rr.numerator());
                    return Value(BigRat(dm.second));
                }
            }
            BigFloat lf = l.to_float();
            BigFloat rf = r.to_float();
            BigFloat q = lf / rf;
            std::string qs = q.to_string(0);
            size_t dot = qs.find('.');
            if (dot != std::string::npos) qs = qs.substr(0, dot);
            BigFloat qi(qs);
            return Value(lf - qi * rf);
        }
        default: return Value::make_error("Unknown op in substitution");
        }
    }
    case ASTNode::UNARYOP: {
        Value operand = substitute(node->left, var, val);
        if (node->op == '-') return -operand;
        return operand;
    }
    case ASTNode::FUNCTION: {
        std::vector<Value> args;
        for (auto& a : node->args) args.push_back(substitute(a, var, val));
        Evaluator temp;
        temp.angle_mode_ = angle_mode_;
        if (node->name == "sin") return temp.eval_sin(args[0]);
        if (node->name == "cos") return temp.eval_cos(args[0]);
        if (node->name == "tan") return temp.eval_tan(args[0]);
        if (node->name == "sqrt") return temp.eval_sqrt(args[0]);
        if (node->name == "abs") return temp.eval_abs(args[0]);
        if (node->name == "ln") return temp.eval_ln(args[0]);
        if (node->name == "lg") return temp.eval_lg(args[0]);
        if (node->name == "exp") return Value(BigFloat::exp(args[0].to_float()));
        return Value(BigFloat::sin_val(args[0].to_float()));
    }
    case ASTNode::FACTORIAL:
        return eval_fact(substitute(node->left, var, val));
    case ASTNode::VEC_LITERAL: {
        std::vector<Value> components;
        for (auto& arg : node->args) {
            Value v = substitute(arg, var, val);
            if (v.is_error()) return v;
            components.push_back(v);
        }
        return Value::make_vector(components);
    }
    case ASTNode::MAT_LITERAL: {
        std::vector<std::vector<Value>> rows;
        for (auto& row : node->mat_rows) {
            std::vector<Value> r;
            for (auto& arg : row) {
                Value v = substitute(arg, var, val);
                if (v.is_error()) return v;
                r.push_back(v);
            }
            rows.push_back(r);
        }
        return Value::make_matrix(rows);
    }
    }
    return Value::make_error("Cannot substitute");
}

Value Evaluator::eval_int(ASTPtr node) {
    if (node->args.size() != 4)
        return Value::make_error("int(expr, var, a, b) requires 4 arguments");
    auto expr = node->args[0];
    std::string var;
    if (node->args[1]->type == ASTNode::VARIABLE) var = node->args[1]->name;
    else return Value::make_error("Second argument of int must be a variable");
    Value a_val = eval_node(node->args[2]);
    Value b_val = eval_node(node->args[3]);
    if (a_val.is_error()) return a_val;
    if (b_val.is_error()) return b_val;
    BigFloat a = a_val.to_float();
    BigFloat b = b_val.to_float();
    int n = 1000;
    BigFloat h = (b - a) / BigFloat(n);
    BigFloat sum(0);
    BigFloat fa = substitute(expr, var, Value(a)).to_float();
    BigFloat fb = substitute(expr, var, Value(b)).to_float();
    sum = sum + fa + fb;
    for (int i = 1; i < n; i++) {
        BigFloat x = a + BigFloat(i) * h;
        Value fx = substitute(expr, var, Value(x));
        if (fx.is_error()) continue;
        BigFloat coeff = (i % 2 == 0) ? BigFloat(2) : BigFloat(4);
        sum = sum + coeff * fx.to_float();
    }
    return Value(sum * h / BigFloat(3));
}

Value Evaluator::eval_diff(ASTPtr node) {
    if (node->args.size() != 3)
        return Value::make_error("diff(expr, var, point) requires 3 arguments");
    auto expr = node->args[0];
    std::string var;
    if (node->args[1]->type == ASTNode::VARIABLE) var = node->args[1]->name;
    else return Value::make_error("Second argument of diff must be a variable");
    Value pt_val = eval_node(node->args[2]);
    if (pt_val.is_error()) return pt_val;
    BigFloat x = pt_val.to_float();
    int p = g_precision;
    BigFloat h("1e-" + std::to_string(p / 3 + 5));
    Value fph = substitute(expr, var, Value(x + h));
    Value fmh = substitute(expr, var, Value(x - h));
    Value f2ph = substitute(expr, var, Value(x + h + h));
    Value f2mh = substitute(expr, var, Value(x - h - h));
    if (fph.is_error()) return fph;
    if (fmh.is_error()) return fmh;
    BigFloat fp = fph.to_float();
    BigFloat fm = fmh.to_float();
    BigFloat f2p = f2ph.to_float();
    BigFloat f2m = f2mh.to_float();
    BigFloat deriv = (BigFloat(-1) * f2p + BigFloat(8) * fp - BigFloat(8) * fm + f2m) / (h * BigFloat(12));
    return Value(deriv);
}

Value Evaluator::eval_sum(ASTPtr node) {
    if (node->args.size() != 4)
        return Value::make_error("sum(expr, var, from, to) requires 4 arguments");
    auto expr = node->args[0];
    std::string var;
    if (node->args[1]->type == ASTNode::VARIABLE) var = node->args[1]->name;
    else return Value::make_error("Second argument of sum must be a variable");
    Value from_val = eval_node(node->args[2]);
    Value to_val = eval_node(node->args[3]);
    if (from_val.is_error()) return from_val;
    if (to_val.is_error()) return to_val;
    if (!from_val.is_rational() || !to_val.is_rational())
        return Value::make_error("Sum bounds must be integers");
    BigRat from_r = from_val.to_rational();
    BigRat to_r = to_val.to_rational();
    if (from_r.denominator() != BigInt(1) || to_r.denominator() != BigInt(1))
        return Value::make_error("Sum bounds must be integers");
    int64_t from = from_r.numerator().to_int64();
    int64_t to = to_r.numerator().to_int64();
    Value result(BigRat(0));
    for (int64_t i = from; i <= to; i++) {
        Value term = substitute(expr, var, Value(BigRat(i)));
        if (term.is_error()) return term;
        result = result + term;
    }
    return result;
}

Value Evaluator::eval_prod(ASTPtr node) {
    if (node->args.size() != 4)
        return Value::make_error("prod(expr, var, from, to) requires 4 arguments");
    auto expr = node->args[0];
    std::string var;
    if (node->args[1]->type == ASTNode::VARIABLE) var = node->args[1]->name;
    else return Value::make_error("Second argument of prod must be a variable");
    Value from_val = eval_node(node->args[2]);
    Value to_val = eval_node(node->args[3]);
    if (from_val.is_error()) return from_val;
    if (to_val.is_error()) return to_val;
    if (!from_val.is_rational() || !to_val.is_rational())
        return Value::make_error("Product bounds must be integers");
    BigRat from_r = from_val.to_rational();
    BigRat to_r = to_val.to_rational();
    if (from_r.denominator() != BigInt(1) || to_r.denominator() != BigInt(1))
        return Value::make_error("Product bounds must be integers");
    int64_t from = from_r.numerator().to_int64();
    int64_t to = to_r.numerator().to_int64();
    Value result(BigRat(1));
    for (int64_t i = from; i <= to; i++) {
        Value term = substitute(expr, var, Value(BigRat(i)));
        if (term.is_error()) return term;
        result = result * term;
    }
    return result;
}

Value Evaluator::eval_complex(const Value& r, const Value& i) {
    return Value::make_complex(r, i);
}

Value Evaluator::eval_re(const Value& v) {
    return v.real_part();
}

Value Evaluator::eval_im(const Value& v) {
    return v.imag_part();
}

Value Evaluator::eval_conj(const Value& v) {
    return v.conjugate();
}

Value Evaluator::eval_arg(const Value& v) {
    return v.argument();
}

Value Evaluator::eval_mod(const Value& v) {
    return v.magnitude();
}

Value Evaluator::eval_rand() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> dist(0, 999999999999999999LL);
    int64_t raw = dist(gen);
    BigFloat numerator(raw);
    BigFloat denominator("1000000000000000000");
    return Value(numerator / denominator);
}

Value Evaluator::eval_randint(const Value& a, const Value& b) {
    if (!a.is_rational() || !b.is_rational())
        return Value::make_error("randint requires integer arguments");
    BigRat ar = a.to_rational(), br = b.to_rational();
    if (ar.denominator() != BigInt(1) || br.denominator() != BigInt(1))
        return Value::make_error("randint requires integer arguments");
    int64_t av = ar.numerator().to_int64();
    int64_t bv = br.numerator().to_int64();
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> dist(av, bv);
    return Value(BigRat(dist(gen)));
}

static std::map<std::string, std::map<std::string, std::string>> build_unit_table() {
    std::map<std::string, std::map<std::string, std::string>> table;
    auto& length = table["length"];
    length["m"] = "1";
    length["km"] = "1000";
    length["cm"] = "0.01";
    length["mm"] = "0.001";
    length["um"] = "0.000001";
    length["nm"] = "0.000000001";
    length["in"] = "0.0254";
    length["ft"] = "0.3048";
    length["yd"] = "0.9144";
    length["mi"] = "1609.344";
    length["nmi"] = "1852";
    length["li"] = "500";
    length["au"] = "149600000000";
    auto& mass = table["mass"];
    mass["kg"] = "1";
    mass["g"] = "0.001";
    mass["mg"] = "0.000001";
    mass["t"] = "1000";
    mass["lb"] = "0.45359237";
    mass["oz"] = "0.028349523125";
    mass["jin"] = "0.5";
    mass["liang"] = "0.05";
    auto& time = table["time"];
    time["s"] = "1";
    time["ms"] = "0.001";
    time["min"] = "60";
    time["h"] = "3600";
    time["d"] = "86400";
    time["wk"] = "604800";
    time["yr"] = "31536000";
    auto& area = table["area"];
    area["m2"] = "1";
    area["km2"] = "1000000";
    area["cm2"] = "0.0001";
    area["ha"] = "10000";
    area["acre"] = "4046.8564224";
    area["mu"] = "666.666666667";
    auto& volume = table["volume"];
    volume["m3"] = "1";
    volume["L"] = "0.001";
    volume["mL"] = "0.000001";
    volume["gal"] = "0.003785411784";
    volume["qt"] = "0.000946352946";
    volume["cup"] = "0.0002365882365";
    auto& temp = table["temperature"];
    temp["C"] = "1";
    temp["K"] = "1";
    temp["F"] = "1";
    auto& speed = table["speed"];
    speed["m/s"] = "1";
    speed["km/h"] = "0.277777777778";
    speed["mph"] = "0.44704";
    speed["kn"] = "0.514444";
    speed["c"] = "299792458";
    auto& data = table["data"];
    data["B"] = "1";
    data["KB"] = "1024";
    data["MB"] = "1048576";
    data["GB"] = "1073741824";
    data["TB"] = "1099511627776";
    data["bit"] = "0.125";
    return table;
}

Value Evaluator::eval_convert(const Value& v, const std::string& from_raw, const std::string& to_raw) {
    static auto table = build_unit_table();
    std::string from = from_raw, to = to_raw;
    auto make_upper = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        return s;
    };
    std::string from_upper = make_upper(from);
    std::string to_upper = make_upper(to);
    if (from_upper == "C" && to_upper == "F") {
        BigFloat val = v.to_float();
        return Value(val * BigFloat(9) / BigFloat(5) + BigFloat(32));
    }
    if (from_upper == "F" && to_upper == "C") {
        BigFloat val = v.to_float();
        return Value((val - BigFloat(32)) * BigFloat(5) / BigFloat(9));
    }
    if (from_upper == "C" && to_upper == "K") {
        return Value(v.to_float() + BigFloat("273.15"));
    }
    if (from_upper == "K" && to_upper == "C") {
        return Value(v.to_float() - BigFloat("273.15"));
    }
    if (from_upper == "F" && to_upper == "K") {
        BigFloat val = v.to_float();
        return Value((val - BigFloat(32)) * BigFloat(5) / BigFloat(9) + BigFloat("273.15"));
    }
    if (from_upper == "K" && to_upper == "F") {
        BigFloat val = v.to_float();
        return Value((val - BigFloat("273.15")) * BigFloat(9) / BigFloat(5) + BigFloat(32));
    }
    for (auto& [category, units] : table) {
        if (category == "temperature") continue;
        auto it_from = units.find(from);
        auto it_to = units.find(to);
        if (it_from != units.end() && it_to != units.end()) {
            BigFloat val = v.to_float();
            BigFloat from_factor(it_from->second);
            BigFloat to_factor(it_to->second);
            BigFloat in_base = val * from_factor;
            return Value(in_base / to_factor);
        }
    }
    return Value::make_error("Unknown unit conversion: " + from + " -> " + to);
}

Value Evaluator::eval_function(ASTPtr node) {
    const std::string& name = node->name;
    auto& args = node->args;

    if (name == "sin") {
        if (args.size() != 1) return Value::make_error("sin requires 1 argument");
        return eval_sin(eval_node(args[0]));
    }
    if (name == "cos") {
        if (args.size() != 1) return Value::make_error("cos requires 1 argument");
        return eval_cos(eval_node(args[0]));
    }
    if (name == "tan") {
        if (args.size() != 1) return Value::make_error("tan requires 1 argument");
        return eval_tan(eval_node(args[0]));
    }
    if (name == "asin") {
        if (args.size() != 1) return Value::make_error("asin requires 1 argument");
        return eval_asin(eval_node(args[0]));
    }
    if (name == "acos") {
        if (args.size() != 1) return Value::make_error("acos requires 1 argument");
        return eval_acos(eval_node(args[0]));
    }
    if (name == "atan") {
        if (args.size() != 1) return Value::make_error("atan requires 1 argument");
        return eval_atan(eval_node(args[0]));
    }
    if (name == "sinh") {
        if (args.size() != 1) return Value::make_error("sinh requires 1 argument");
        return eval_sinh(eval_node(args[0]));
    }
    if (name == "cosh") {
        if (args.size() != 1) return Value::make_error("cosh requires 1 argument");
        return eval_cosh(eval_node(args[0]));
    }
    if (name == "tanh") {
        if (args.size() != 1) return Value::make_error("tanh requires 1 argument");
        return eval_tanh(eval_node(args[0]));
    }
    if (name == "sqrt") {
        if (args.size() != 1) return Value::make_error("sqrt requires 1 argument");
        return eval_sqrt(eval_node(args[0]));
    }
    if (name == "cbrt") {
        if (args.size() != 1) return Value::make_error("cbrt requires 1 argument");
        return eval_cbrt(eval_node(args[0]));
    }
    if (name == "nrt") {
        if (args.size() != 2) return Value::make_error("nrt(n, x) requires 2 arguments");
        return eval_nrt(eval_node(args[1]), eval_node(args[0]));
    }
    if (name == "abs") {
        if (args.size() != 1) return Value::make_error("abs requires 1 argument");
        return eval_abs(eval_node(args[0]));
    }
    if (name == "ln") {
        if (args.size() != 1) return Value::make_error("ln requires 1 argument");
        return eval_ln(eval_node(args[0]));
    }
    if (name == "lg") {
        if (args.size() != 1) return Value::make_error("lg requires 1 argument");
        return eval_lg(eval_node(args[0]));
    }
    if (name == "log") {
        if (args.size() == 1) return eval_ln(eval_node(args[0]));
        if (args.size() == 2) return eval_log(eval_node(args[0]), eval_node(args[1]));
        return Value::make_error("log requires 1 or 2 arguments");
    }
    if (name == "exp") {
        if (args.size() != 1) return Value::make_error("exp requires 1 argument");
        Value arg = eval_node(args[0]);
        if (arg.type == Value::SURDS && arg.surds.is_rational()) {
            BigRat r = arg.surds.to_rational();
            if (r == BigRat(1)) return Value(SurdsExpr(BigRat(1), BigInt(-2)));
            if (r == BigRat(0)) return Value(BigRat(1));
        }
        return Value(BigFloat::exp(arg.to_float()));
    }
    if (name == "fact") {
        if (args.size() != 1) return Value::make_error("fact requires 1 argument");
        return eval_fact(eval_node(args[0]));
    }
    if (name == "perm" || name == "P") {
        if (args.size() != 2) return Value::make_error("perm(n, k) requires 2 arguments");
        return eval_perm(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "comb" || name == "C") {
        if (args.size() != 2) return Value::make_error("comb(n, k) requires 2 arguments");
        return eval_comb(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "gcd") {
        if (args.size() != 2) return Value::make_error("gcd requires 2 arguments");
        return eval_gcd(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "lcm") {
        if (args.size() != 2) return Value::make_error("lcm requires 2 arguments");
        return eval_lcm(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "factor") {
        if (args.size() != 1) return Value::make_error("factor requires 1 argument");
        return eval_factor(eval_node(args[0]));
    }
    if (name == "int") {
        return eval_int(node);
    }
    if (name == "diff") {
        return eval_diff(node);
    }
    if (name == "sum") {
        return eval_sum(node);
    }
    if (name == "prod") {
        return eval_prod(node);
    }
    if (name == "complex") {
        if (args.size() != 2) return Value::make_error("complex(a, b) requires 2 arguments");
        return eval_complex(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "re") {
        if (args.size() != 1) return Value::make_error("re requires 1 argument");
        return eval_re(eval_node(args[0]));
    }
    if (name == "im") {
        if (args.size() != 1) return Value::make_error("im requires 1 argument");
        return eval_im(eval_node(args[0]));
    }
    if (name == "conj") {
        if (args.size() != 1) return Value::make_error("conj requires 1 argument");
        return eval_conj(eval_node(args[0]));
    }
    if (name == "arg") {
        if (args.size() != 1) return Value::make_error("arg requires 1 argument");
        return eval_arg(eval_node(args[0]));
    }
    if (name == "mod") {
        if (args.size() != 1) return Value::make_error("mod requires 1 argument");
        return eval_mod(eval_node(args[0]));
    }
    if (name == "rand") {
        return eval_rand();
    }
    if (name == "randint") {
        if (args.size() != 2) return Value::make_error("randint(a, b) requires 2 arguments");
        return eval_randint(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "convert") {
        if (args.size() != 3) return Value::make_error("convert(val, from, to) requires 3 arguments");
        Value val = eval_node(args[0]);
        std::string from, to;
        if (args[1]->type == ASTNode::VARIABLE || args[1]->type == ASTNode::CONSTANT)
            from = args[1]->name;
        else if (args[1]->type == ASTNode::NUMBER)
            from = args[1]->number.to_string();
        else return Value::make_error("Invalid 'from' unit");
        if (args[2]->type == ASTNode::VARIABLE || args[2]->type == ASTNode::CONSTANT)
            to = args[2]->name;
        else if (args[2]->type == ASTNode::NUMBER)
            to = args[2]->number.to_string();
        else return Value::make_error("Invalid 'to' unit");
        return eval_convert(val, from, to);
    }
    if (name == "solve") {
        return eval_solve(node->args);
    }
    if (name == "deg") {
        if (args.size() == 1) {
            Evaluator temp;
            temp.angle_mode_ = DEG;
            return temp.eval_node(args[0]);
        }
        angle_mode_ = DEG;
        return Value(BigRat(0));
    }
    if (name == "rad") {
        if (args.size() == 1) {
            Evaluator temp;
            temp.angle_mode_ = RAD;
            return temp.eval_node(args[0]);
        }
        angle_mode_ = RAD;
        return Value(BigRat(0));
    }
    if (name == "floor") {
        if (args.size() != 1) return Value::make_error("floor requires 1 argument");
        Value v = eval_node(args[0]);
        if (v.is_rational()) {
            BigRat r = v.to_rational();
            BigInt q = r.numerator() / r.denominator();
            if (r < BigRat(0) && r.numerator() % r.denominator() != BigInt(0))
                q = q - BigInt(1);
            return Value(q);
        }
        return Value(v.to_float());
    }
    if (name == "ceil") {
        if (args.size() != 1) return Value::make_error("ceil requires 1 argument");
        Value v = eval_node(args[0]);
        if (v.is_rational()) {
            BigRat r = v.to_rational();
            BigInt q = r.numerator() / r.denominator();
            if (r > BigRat(0) && r.numerator() % r.denominator() != BigInt(0))
                q = q + BigInt(1);
            return Value(q);
        }
        return Value(v.to_float());
    }
    if (name == "round") {
        if (args.size() != 1) return Value::make_error("round requires 1 argument");
        Value v = eval_node(args[0]);
        if (v.is_rational()) {
            BigRat r = v.to_rational();
            BigRat half(1, 2);
            if (r >= BigRat(0)) {
                BigRat rounded = r + half;
                return Value(rounded.numerator() / rounded.denominator());
            } else {
                BigRat rounded = r - half;
                BigInt q = rounded.numerator() / rounded.denominator();
                return Value(q);
            }
        }
        BigFloat f = v.to_float();
        return Value(f.round());
    }
    if (name == "sign") {
        if (args.size() != 1) return Value::make_error("sign requires 1 argument");
        Value v = eval_node(args[0]);
        if (v.is_rational()) return Value(BigRat(v.to_rational().sign()));
        return Value(BigRat(v.to_float().sign()));
    }
    if (name == "max") {
        if (args.size() < 2) return Value::make_error("max requires at least 2 arguments");
        Value result = eval_node(args[0]);
        for (size_t i = 1; i < args.size(); i++) {
            Value v = eval_node(args[i]);
            if (v.to_float() > result.to_float()) result = v;
        }
        return result;
    }
    if (name == "min") {
        if (args.size() < 2) return Value::make_error("min requires at least 2 arguments");
        Value result = eval_node(args[0]);
        for (size_t i = 1; i < args.size(); i++) {
            Value v = eval_node(args[i]);
            if (v.to_float() < result.to_float()) result = v;
        }
        return result;
    }
    if (name == "pow") {
        if (args.size() != 2) return Value::make_error("pow requires 2 arguments");
        Value base = eval_node(args[0]);
        Value exp = eval_node(args[1]);
        return Value(BigFloat::pow_val(base.to_float(), exp.to_float()));
    }
    if (name == "root") {
        if (args.size() != 2) return Value::make_error("root requires 2 arguments");
        Value n = eval_node(args[0]);
        Value x = eval_node(args[1]);
        return eval_nrt(x, n);
    }
    if (name == "log2") {
        if (args.size() != 1) return Value::make_error("log2 requires 1 argument");
        return Value(BigFloat::log2(eval_node(args[0]).to_float()));
    }
    if (name == "atan2") {
        if (args.size() != 2) return Value::make_error("atan2 requires 2 arguments");
        Value y = eval_node(args[0]);
        Value x = eval_node(args[1]);
        return Value(BigFloat::atan_val(y.to_float() / x.to_float()));
    }
    if (name == "hypot") {
        if (args.size() != 2) return Value::make_error("hypot requires 2 arguments");
        Value a = eval_node(args[0]);
        Value b = eval_node(args[1]);
        return Value(BigFloat::sqrt_val(
            a.to_float() * a.to_float() + b.to_float() * b.to_float()));
    }
    if (name == "vecmod") {
        if (args.size() != 1) return Value::make_error("vecmod requires 1 argument");
        return eval_vecmod(eval_node(args[0]));
    }
    if (name == "dot") {
        if (args.size() != 2) return Value::make_error("dot requires 2 arguments");
        return eval_dot(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "cross") {
        if (args.size() != 2) return Value::make_error("cross requires 2 arguments");
        return eval_cross(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "scalarmul") {
        if (args.size() != 2) return Value::make_error("scalarmul requires 2 arguments");
        return eval_scalarmul(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "mixed") {
        if (args.size() != 3) return Value::make_error("mixed requires 3 arguments");
        return eval_mixed(eval_node(args[0]), eval_node(args[1]), eval_node(args[2]));
    }
    if (name == "proj") {
        if (args.size() != 2) return Value::make_error("proj requires 2 arguments");
        return eval_proj(eval_node(args[0]), eval_node(args[1]));
    }
    if (name == "decompose") {
        if (args.size() == 3)
            return eval_decompose(eval_node(args[0]), eval_node(args[1]), eval_node(args[2]));
        if (args.size() == 4)
            return eval_decompose3d(eval_node(args[0]), eval_node(args[1]), eval_node(args[2]), eval_node(args[3]));
        return Value::make_error("decompose requires 3 or 4 arguments");
    }
    if (name == "solve") {
        return eval_solve(args);
    }
    if (name == "matrix") {
        return eval_matrix(args);
    }
    if (name == "det") {
        if (args.size() != 1) return Value::make_error("det requires 1 argument");
        return eval_det(eval_node(args[0]));
    }
    if (name == "inv") {
        if (args.size() != 1) return Value::make_error("inv requires 1 argument");
        return eval_inv(eval_node(args[0]));
    }
    if (name == "eigen") {
        if (args.size() != 1) return Value::make_error("eigen requires 1 argument");
        return eval_eigen(eval_node(args[0]));
    }
    if (name == "trace") {
        if (args.size() != 1) return Value::make_error("trace requires 1 argument");
        return eval_trace(eval_node(args[0]));
    }
    if (name == "transpose") {
        if (args.size() != 1) return Value::make_error("transpose requires 1 argument");
        return eval_transpose(eval_node(args[0]));
    }
    if (name == "identity" || name == "eye") {
        if (args.size() != 1) return Value::make_error("identity requires 1 argument");
        return eval_identity(eval_node(args[0]));
    }
    if (name == "mean") {
        std::vector<Value> vals;
        for (auto& a : args) { Value v = eval_node(a); if (v.is_error()) return v; vals.push_back(v); }
        return eval_mean(vals);
    }
    if (name == "stddev") {
        std::vector<Value> vals;
        for (auto& a : args) { Value v = eval_node(a); if (v.is_error()) return v; vals.push_back(v); }
        return eval_stddev(vals);
    }
    if (name == "variance") {
        std::vector<Value> vals;
        for (auto& a : args) { Value v = eval_node(a); if (v.is_error()) return v; vals.push_back(v); }
        return eval_variance(vals);
    }
    if (name == "median") {
        std::vector<Value> vals;
        for (auto& a : args) { Value v = eval_node(a); if (v.is_error()) return v; vals.push_back(v); }
        return eval_median(vals);
    }
    if (name == "simplify") {
        if (args.size() != 1) return Value::make_error("simplify requires 1 argument");
        return eval_simplify(args[0]);
    }
    return Value::make_error("Unknown function: " + name);
}

Value Evaluator::eval_vecmod(const Value& v) {
    if (!v.is_vector())
        return Value::make_error("vecmod requires a vector");
    if (v.vec_dim() == 0)
        return Value::make_error("vecmod: empty vector");
    Value sum(BigRat(0));
    for (auto& c : v.vec) {
        sum = sum + c * c;
    }
    return eval_sqrt(sum);
}

Value Evaluator::eval_dot(const Value& a, const Value& b) {
    if (!a.is_vector() || !b.is_vector())
        return Value::make_error("dot requires two vectors");
    if (a.vec_dim() != b.vec_dim())
        return Value::make_error("dot: vectors must have same dimension");
    if (a.vec_dim() == 0)
        return Value::make_error("dot: empty vectors");
    Value sum(BigRat(0));
    for (int i = 0; i < a.vec_dim(); i++) {
        sum = sum + a.vec[i] * b.vec[i];
    }
    return sum;
}

Value Evaluator::eval_cross(const Value& a, const Value& b) {
    if (!a.is_vector() || !b.is_vector())
        return Value::make_error("cross requires two vectors");
    if (a.vec_dim() != 3 || b.vec_dim() != 3)
        return Value::make_error("cross product requires 3D vectors");
    Value c1 = a.vec[1] * b.vec[2] - a.vec[2] * b.vec[1];
    Value c2 = a.vec[2] * b.vec[0] - a.vec[0] * b.vec[2];
    Value c3 = a.vec[0] * b.vec[1] - a.vec[1] * b.vec[0];
    return Value::make_vector({c1, c2, c3});
}

Value Evaluator::eval_scalarmul(const Value& s, const Value& v) {
    if (!v.is_vector())
        return Value::make_error("scalarmul requires a scalar and a vector");
    if (v.vec_dim() == 0)
        return Value::make_error("scalarmul: empty vector");
    return s * v;
}

Value Evaluator::eval_mixed(const Value& a, const Value& b, const Value& c) {
    if (!a.is_vector() || !b.is_vector() || !c.is_vector())
        return Value::make_error("mixed requires three vectors");
    if (a.vec_dim() != 3 || b.vec_dim() != 3 || c.vec_dim() != 3)
        return Value::make_error("mixed product requires 3D vectors");
    Value cp = eval_cross(b, c);
    if (cp.is_error()) return cp;
    return eval_dot(a, cp);
}

Value Evaluator::eval_proj(const Value& a, const Value& b) {
    if (!a.is_vector() || !b.is_vector())
        return Value::make_error("proj requires two vectors");
    if (a.vec_dim() != b.vec_dim())
        return Value::make_error("proj: vectors must have same dimension");
    if (b.is_zero())
        return Value::make_error("proj: cannot project onto zero vector");
    if (a.vec_dim() == 0)
        return Value::make_error("proj: empty vectors");
    Value d = eval_dot(a, b);
    if (d.is_error()) return d;
    Value bb = eval_dot(b, b);
    if (bb.is_error()) return bb;
    Value coeff = d / bb;
    return eval_scalarmul(coeff, b);
}

Value Evaluator::eval_decompose(const Value& a, const Value& b, const Value& c) {
    if (!a.is_vector() || !b.is_vector() || !c.is_vector())
        return Value::make_error("decompose requires vectors");
    if (a.vec_dim() != b.vec_dim() || b.vec_dim() != c.vec_dim())
        return Value::make_error("decompose: vectors must have same dimension");
    if (b.is_zero() || c.is_zero())
        return Value::make_error("decompose: basis vectors cannot be zero");

    if (a.vec_dim() == 2) {
        Value det = b.vec[0] * c.vec[1] - b.vec[1] * c.vec[0];
        if (det.is_zero())
            return Value::make_error("decompose: basis vectors are collinear");
        Value alpha = (a.vec[0] * c.vec[1] - a.vec[1] * c.vec[0]) / det;
        Value beta = (b.vec[0] * a.vec[1] - b.vec[1] * a.vec[0]) / det;
        return Value::make_vector({alpha, beta});
    }

    if (a.vec_dim() == 3)
        return Value::make_error("decompose: for 3D, use decompose(a, b, c, d) with 3 basis vectors");

    return Value::make_error("decompose: only 2D and 3D vectors supported");
}

Value Evaluator::eval_decompose3d(const Value& a, const Value& b, const Value& c, const Value& d) {
    if (!a.is_vector() || !b.is_vector() || !c.is_vector() || !d.is_vector())
        return Value::make_error("decompose requires vectors");
    if (a.vec_dim() != b.vec_dim() || b.vec_dim() != c.vec_dim() || c.vec_dim() != d.vec_dim())
        return Value::make_error("decompose: vectors must have same dimension");
    if (a.vec_dim() != 3)
        return Value::make_error("decompose3d requires 3D vectors");

    Value det = b.vec[0] * (c.vec[1] * d.vec[2] - c.vec[2] * d.vec[1])
              - c.vec[0] * (b.vec[1] * d.vec[2] - b.vec[2] * d.vec[1])
              + d.vec[0] * (b.vec[1] * c.vec[2] - b.vec[2] * c.vec[1]);
    if (det.is_zero())
        return Value::make_error("decompose: basis vectors are coplanar");

    Value da = a.vec[0] * (c.vec[1] * d.vec[2] - c.vec[2] * d.vec[1])
             - c.vec[0] * (a.vec[1] * d.vec[2] - a.vec[2] * d.vec[1])
             + d.vec[0] * (a.vec[1] * c.vec[2] - a.vec[2] * c.vec[1]);

    Value db = b.vec[0] * (a.vec[1] * d.vec[2] - a.vec[2] * d.vec[1])
             - a.vec[0] * (b.vec[1] * d.vec[2] - b.vec[2] * d.vec[1])
             + d.vec[0] * (b.vec[1] * a.vec[2] - b.vec[2] * a.vec[1]);

    Value dc = b.vec[0] * (c.vec[1] * a.vec[2] - c.vec[2] * a.vec[1])
             - c.vec[0] * (b.vec[1] * a.vec[2] - b.vec[2] * a.vec[1])
             + a.vec[0] * (b.vec[1] * c.vec[2] - b.vec[2] * c.vec[1]);

    Value alpha = da / det;
    Value beta = db / det;
    Value gamma = dc / det;
    return Value::make_vector({alpha, beta, gamma});
}

Value Evaluator::eval_solve_quadratic(const Value& a, const Value& b, const Value& c) {
    if (a.is_zero()) {
        if (b.is_zero()) {
            if (c.is_zero()) return Value::make_error("Infinite solutions (0=0)");
            return Value::make_error("No solution (constant = 0)");
        }
        Value root = -c / b;
        return Value::make_string("Root: " + format_result(root));
    }
    Value D = b * b - Value(BigRat(4)) * a * c;
    Value sqrtD = eval_sqrt(D);
    if (sqrtD.is_error()) return sqrtD;
    Value two_a = Value(BigRat(2)) * a;
    Value x1 = (-b + sqrtD) / two_a;
    Value x2 = (-b - sqrtD) / two_a;

    Value vx = -b / (Value(BigRat(2)) * a);
    Value vy = a * vx * vx + b * vx + c;

    std::string result;
    result += "Roots: " + format_result(x1);
    if (format_result(x1) != format_result(x2))
        result += ", " + format_result(x2);
    result += "\n";
    result += "Delta: " + format_result(D) + "\n";

    Value a_float = a;
    if (a_float.type == Value::SURDS) a_float = Value(a_float.to_float());
    bool a_positive = std::stod(a_float.to_string()) > 0;

    if (a_positive) {
        result += "Min: " + format_result(vy) + " at x = " + format_result(vx);
    } else {
        result += "Max: " + format_result(vy) + " at x = " + format_result(vx);
    }

    return Value::make_string(result);
}

static double poly_eval(const std::vector<double>& coeffs, double x) {
    double result = 0;
    for (size_t i = 0; i < coeffs.size(); i++) {
        result = result * x + coeffs[i];
    }
    return result;
}

static std::vector<double> poly_deriv(const std::vector<double>& coeffs) {
    int n = (int)coeffs.size() - 1;
    std::vector<double> d;
    for (int i = 0; i < n; i++) {
        d.push_back(coeffs[i] * (n - i));
    }
    return d;
}

static std::string fmt_d(double v) {
    if (fabs(v) < 5e-10) return "0";
    if (fabs(v - round(v)) < 1e-9) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(0) << round(v);
        return oss.str();
    }
    std::ostringstream oss;
    oss << std::setprecision(10) << v;
    std::string s = oss.str();
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    return s;
}

static std::string fmt_complex(const std::complex<double>& z) {
    if (fabs(z.imag()) < 1e-10) return fmt_d(z.real());
    if (fabs(z.real()) < 1e-10) {
        if (fabs(z.imag() - 1.0) < 1e-10) return "i";
        if (fabs(z.imag() + 1.0) < 1e-10) return "-i";
        return fmt_d(z.imag()) + "*i";
    }
    std::string s = fmt_d(z.real()) + " + " + fmt_d(z.imag()) + "*i";
    return s;
}

static std::string format_extrema(const std::vector<double>& coeffs,
                                   const std::vector<std::complex<double>>& crit_pts) {
    struct Extremum {
        double x, y;
        bool is_max;
    };
    std::vector<Extremum> extrema;
    auto d1 = poly_deriv(coeffs);
    auto d2 = poly_deriv(d1);
    for (auto& cp : crit_pts) {
        if (fabs(cp.imag()) > 1e-8) continue;
        double x = cp.real();
        double fpp = poly_eval(d2, x);
        if (fabs(fpp) < 1e-12) continue;
        double y = poly_eval(coeffs, x);
        extrema.push_back({x, y, fpp < 0});
    }
    if (extrema.empty()) return "No local extrema";
    std::sort(extrema.begin(), extrema.end(),
              [](const Extremum& a, const Extremum& b) { return a.y < b.y; });
    std::string result;
    for (auto& e : extrema) {
        if (!result.empty()) result += "\n";
        if (e.is_max)
            result += "Local max: " + fmt_d(e.y) + " at x = " + fmt_d(e.x);
        else
            result += "Local min: " + fmt_d(e.y) + " at x = " + fmt_d(e.x);
    }
    int deg = (int)coeffs.size() - 1;
    if (deg == 4) {
        bool a_pos = coeffs[0] > 0;
        if (a_pos) {
            auto& gmin = extrema.front();
            result += "\nGlobal min: " + fmt_d(gmin.y) + " at x = " + fmt_d(gmin.x);
        } else {
            auto& gmax = extrema.back();
            result += "\nGlobal max: " + fmt_d(gmax.y) + " at x = " + fmt_d(gmax.x);
        }
    }
    return result;
}

Value Evaluator::eval_solve(const std::vector<ASTPtr>& args) {
    int n = (int)args.size();
    if (n < 3 || n > 5)
        return Value::make_error("solve requires 3-5 arguments (quadratic to quartic)");

    if (n == 3) {
        Value a = eval_node(args[0]);
        Value b = eval_node(args[1]);
        Value c = eval_node(args[2]);
        if (a.is_error()) return a;
        if (b.is_error()) return b;
        if (c.is_error()) return c;
        return eval_solve_quadratic(a, b, c);
    }

    std::vector<double> coeffs;
    for (auto& arg : args) {
        Value v = eval_node(arg);
        if (v.is_error()) return v;
        coeffs.push_back(std::stod(v.to_string()));
    }

    if (n == 4) {
        auto roots = solve_cubic_d(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
        std::string result = "Roots: ";
        bool first = true;
        for (auto& r : roots) {
            if (!first) result += ", ";
            first = false;
            result += fmt_complex(r);
        }
        double da = 3 * coeffs[0], db = 2 * coeffs[1], dc = coeffs[2];
        auto crit = solve_quad_d(da, db, dc);
        result += "\n" + format_extrema(coeffs, crit);
        return Value::make_string(result);
    }

    auto roots = solve_quartic_d(coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]);
    std::string result = "Roots: ";
    bool first = true;
    for (auto& r : roots) {
        if (!first) result += ", ";
        first = false;
        result += fmt_complex(r);
    }
    double da = 4 * coeffs[0], db = 3 * coeffs[1], dc = 2 * coeffs[2], dd = coeffs[3];
    auto crit = solve_cubic_d(da, db, dc, dd);
    result += "\n" + format_extrema(coeffs, crit);
    return Value::make_string(result);
}

std::string Evaluator::format_result(const Value& v, int base) {
    if (v.is_error()) return v.to_string();
    if (v.is_string()) return v.to_string();
    if (v.is_vector()) {
        std::string r = "(";
        for (size_t i = 0; i < v.vec.size(); i++) {
            if (i > 0) r += ", ";
            r += format_result(v.vec[i], base);
        }
        r += ")";
        return r;
    }
    if (v.type == Value::SURDS) {
        if (v.surds.is_rational()) {
            BigRat r = v.surds.to_rational();
            if (r.denominator() == BigInt(1)) {
                if (base != 10) {
                    return r.numerator().to_base_string(base);
                }
                return r.numerator().to_string();
            }
            if (base != 10) {
                BigFloat fv = r.to_decimal_string().empty() ? v.to_float() : BigFloat(r);
                return fv.to_string();
            }
            std::string frac = r.to_string();
            std::string dec = r.to_decimal_string();
            return frac + " (= " + dec + ")";
        }
        std::string exact = v.surds.to_string();
        std::string approx = v.to_float().to_string();
        return exact + " ~= " + approx;
    }
    return v.to_string();
}
