#ifndef CCALC_H
#define CCALC_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cassert>
#include <functional>
#include <random>
#include <optional>
#include <variant>
#include <stdexcept>
#include <cctype>
#include <utility>
#include <numeric>
#include <iomanip>

extern int g_precision;

class CalcError : public std::runtime_error {
public:
    explicit CalcError(const std::string& msg) : std::runtime_error(msg) {}
};

class BigInt {
public:
    static const int64_t BASE = 1000000000LL;
    static const int BASE_DIGITS = 9;

    BigInt() : neg_(false) {}
    BigInt(int64_t n);
    BigInt(const std::string& s);

    BigInt operator+() const { return *this; }
    BigInt operator-() const;
    BigInt operator+(const BigInt& o) const;
    BigInt operator-(const BigInt& o) const;
    BigInt operator*(const BigInt& o) const;
    BigInt operator/(const BigInt& o) const;
    BigInt operator%(const BigInt& o) const;

    BigInt& operator+=(const BigInt& o) { return *this = *this + o; }
    BigInt& operator-=(const BigInt& o) { return *this = *this - o; }
    BigInt& operator*=(const BigInt& o) { return *this = *this * o; }
    BigInt& operator/=(const BigInt& o) { return *this = *this / o; }
    BigInt& operator%=(const BigInt& o) { return *this = *this % o; }

    bool operator==(const BigInt& o) const;
    bool operator!=(const BigInt& o) const { return !(*this == o); }
    bool operator<(const BigInt& o) const;
    bool operator<=(const BigInt& o) const { return !(o < *this); }
    bool operator>(const BigInt& o) const { return o < *this; }
    bool operator>=(const BigInt& o) const { return !(*this < o); }

    BigInt abs() const;
    int sign() const;
    bool is_zero() const;

    std::string to_string() const;
    std::string to_decimal_string() const;
    std::string to_base_string(int base) const;

    static BigInt from_base_string(const std::string& s, int base);

    static BigInt gcd(BigInt a, BigInt b);
    static BigInt lcm(const BigInt& a, const BigInt& b);
    static BigInt pow(BigInt base, BigInt exp);
    static BigInt factorial(int64_t n);

    friend std::ostream& operator<<(std::ostream& os, const BigInt& n) {
        return os << n.to_string();
    }

    int64_t to_int64() const;

    const std::vector<int64_t>& digits() const { return digits_; }
    bool negative() const { return neg_; }

private:
    std::vector<int64_t> digits_;
    bool neg_;

    void trim();
    BigInt raw_add(const BigInt& o) const;
    BigInt raw_sub(const BigInt& o) const;
    int raw_compare(const BigInt& o) const;

    static std::vector<int64_t> karatsuba(const std::vector<int64_t>& x,
                                           const std::vector<int64_t>& y);
    static std::vector<int64_t> add_vec(const std::vector<int64_t>& a,
                                         const std::vector<int64_t>& b);
    static std::vector<int64_t> sub_vec(const std::vector<int64_t>& a,
                                         const std::vector<int64_t>& b);
    static std::vector<int64_t> shift(const std::vector<int64_t>& v, size_t offset);

public:
    std::pair<BigInt, BigInt> divmod(const BigInt& o) const;
};

class BigRat {
public:
    BigRat() : num_(0), den_(1) {}
    BigRat(int64_t n) : num_(n), den_(1) {}
    BigRat(int64_t n, int64_t d);
    BigRat(const BigInt& n) : num_(n), den_(1) {}
    BigRat(const BigInt& n, const BigInt& d);
    BigRat(const std::string& s);

    BigRat operator+(const BigRat& o) const;
    BigRat operator-(const BigRat& o) const;
    BigRat operator*(const BigRat& o) const;
    BigRat operator/(const BigRat& o) const;
    BigRat operator-() const;

    BigRat& operator+=(const BigRat& o) { return *this = *this + o; }
    BigRat& operator-=(const BigRat& o) { return *this = *this - o; }
    BigRat& operator*=(const BigRat& o) { return *this = *this * o; }
    BigRat& operator/=(const BigRat& o) { return *this = *this / o; }

    bool operator==(const BigRat& o) const;
    bool operator!=(const BigRat& o) const { return !(*this == o); }
    bool operator<(const BigRat& o) const;
    bool operator<=(const BigRat& o) const { return !(o < *this); }
    bool operator>(const BigRat& o) const { return o < *this; }
    bool operator>=(const BigRat& o) const { return !(*this < o); }

    BigRat abs() const;
    int sign() const;
    bool is_zero() const;

    std::string to_string() const;
    std::string to_mixed_string() const;
    std::string to_decimal_string(int precision = 0) const;

    const BigInt& numerator() const { return num_; }
    const BigInt& denominator() const { return den_; }

    static BigRat pow(const BigRat& base, int64_t exp);

    friend std::ostream& operator<<(std::ostream& os, const BigRat& r) {
        return os << r.to_string();
    }

private:
    BigInt num_, den_;
    void simplify();
};

class BigFloat {
public:
    BigFloat() : mantissa_(0), exp_(0) {}
    BigFloat(int n) : BigFloat((int64_t)n) {}
    BigFloat(int64_t n);
    explicit BigFloat(double n);
    explicit BigFloat(const BigInt& n);
    explicit BigFloat(const BigRat& r);
    BigFloat(const std::string& s);

    BigFloat operator+(const BigFloat& o) const;
    BigFloat operator-(const BigFloat& o) const;
    BigFloat operator*(const BigFloat& o) const;
    BigFloat operator/(const BigFloat& o) const;
    BigFloat operator-() const;

    bool operator==(const BigFloat& o) const;
    bool operator!=(const BigFloat& o) const { return !(*this == o); }
    bool operator<(const BigFloat& o) const;
    bool operator<=(const BigFloat& o) const { return !(o < *this); }
    bool operator>(const BigFloat& o) const { return o < *this; }
    bool operator>=(const BigFloat& o) const { return !(*this < o); }

    BigFloat abs() const;
    BigFloat round() const;
    int sign() const;
    bool is_zero() const;

    std::string to_string(int prec = 0) const;

    static BigFloat pi(int prec = 0);
    static BigFloat e_val(int prec = 0);
    static BigFloat ln2(int prec = 0);
    static BigFloat exp(const BigFloat& x, int prec = 0);
    static BigFloat ln(const BigFloat& x, int prec = 0);
    static BigFloat log10(const BigFloat& x, int prec = 0);
    static BigFloat log2(const BigFloat& x, int prec = 0);
    static BigFloat log_base(const BigFloat& x, const BigFloat& base, int prec = 0);
    static BigFloat sin_val(const BigFloat& x, int prec = 0);
    static BigFloat cos_val(const BigFloat& x, int prec = 0);
    static BigFloat tan_val(const BigFloat& x, int prec = 0);
    static BigFloat asin_val(const BigFloat& x, int prec = 0);
    static BigFloat acos_val(const BigFloat& x, int prec = 0);
    static BigFloat atan_val(const BigFloat& x, int prec = 0);
    static BigFloat sinh_val(const BigFloat& x, int prec = 0);
    static BigFloat cosh_val(const BigFloat& x, int prec = 0);
    static BigFloat tanh_val(const BigFloat& x, int prec = 0);
    static BigFloat sqrt_val(const BigFloat& x, int prec = 0);
    static BigFloat cbrt_val(const BigFloat& x, int prec = 0);
    static BigFloat nrt_val(const BigFloat& x, int64_t n, int prec = 0);
    static BigFloat pow_val(const BigFloat& base, const BigFloat& exp, int prec = 0);
    static BigFloat factorial(int64_t n);

    friend std::ostream& operator<<(std::ostream& os, const BigFloat& f) {
        return os << f.to_string();
    }

    const BigInt& mantissa() const { return mantissa_; }
    int64_t exponent() const { return exp_; }

private:
    BigInt mantissa_;
    int64_t exp_;
    static BigFloat pi_cache_;

    void normalize(int prec = 0);
    BigFloat mul_int(const BigInt& n) const;
};

struct SurdsTerm {
    BigRat coeff;
    BigInt radicand;

    SurdsTerm() : radicand(1) {}
    SurdsTerm(const BigRat& c, const BigInt& r) : coeff(c), radicand(r) {}
};

class SurdsExpr {
public:
    std::vector<SurdsTerm> terms;

    SurdsExpr();
    explicit SurdsExpr(const BigRat& r);
    SurdsExpr(const BigRat& c, const BigInt& radicand);

    SurdsExpr operator+(const SurdsExpr& o) const;
    SurdsExpr operator-(const SurdsExpr& o) const;
    SurdsExpr operator*(const SurdsExpr& o) const;
    SurdsExpr operator/(const SurdsExpr& o) const;
    SurdsExpr operator-() const;

    bool operator==(const SurdsExpr& o) const;
    bool is_zero() const;
    bool is_rational() const;
    bool is_positive() const;
    bool is_negative() const;

    BigRat to_rational() const;
    BigFloat to_float(int prec = 0) const;
    std::string to_string() const;

    void simplify();
    void combine_like_terms();

    static BigInt make_square_free(const BigInt& n);
};

class Value;

struct ComplexVal {
    std::shared_ptr<Value> real;
    std::shared_ptr<Value> imag;
    ComplexVal() = default;
    ComplexVal(const Value& r, const Value& i);
};

class Value {
public:
    enum Type { SURDS, FLOAT, COMPLEX, VECTOR, MATRIX, STRING, ERROR };

    Type type;

    SurdsExpr surds;
    BigFloat float_val;
    ComplexVal complex;
    std::vector<Value> vec;
    std::vector<std::vector<Value>> mat;
    int mat_rows() const { return is_matrix() ? (int)mat.size() : 0; }
    int mat_cols() const { return is_matrix() && !mat.empty() ? (int)mat[0].size() : 0; }
    std::string error_msg;

    Value() : type(SURDS), surds(BigRat(0)) {}
    Value(const BigRat& r) : type(SURDS), surds(r) {}
    Value(const BigInt& n) : type(SURDS), surds(BigRat(n)) {}
    Value(int64_t n) : type(SURDS), surds(BigRat(n)) {}
    Value(const SurdsExpr& s) : type(SURDS), surds(s) {}
    Value(const BigFloat& f) : type(FLOAT), float_val(f) {}
    Value(const ComplexVal& c) : type(COMPLEX), complex(c) {}
    Value(Type, const std::string& err) : type(ERROR), error_msg(err) {}

    static Value make_error(const std::string& msg) {
        return Value(ERROR, msg);
    }
    static Value make_string(const std::string& s) {
        Value v;
        v.type = STRING;
        v.error_msg = s;
        return v;
    }
    static Value make_complex(const Value& r, const Value& i);
    static Value make_vector(const std::vector<Value>& components) {
        Value v;
        v.type = VECTOR;
        v.vec = components;
        return v;
    }
    static Value make_matrix(const std::vector<std::vector<Value>>& rows) {
        Value v;
        v.type = MATRIX;
        v.mat = rows;
        return v;
    }

    bool is_error() const { return type == ERROR; }
    bool is_string() const { return type == STRING; }
    bool is_exact() const;
    bool is_rational() const;
    bool is_complex() const { return type == COMPLEX; }
    bool is_vector() const { return type == VECTOR; }
    bool is_matrix() const { return type == MATRIX; }
    bool is_zero() const;

    int vec_dim() const { return is_vector() ? (int)vec.size() : 0; }

    BigRat to_rational() const;
    BigFloat to_float(int prec = 0) const;
    std::string to_string() const;
    std::string to_latex() const;

    Value operator+(const Value& o) const;
    Value operator-(const Value& o) const;
    Value operator*(const Value& o) const;
    Value operator/(const Value& o) const;
    Value operator-() const;

    Value real_part() const;
    Value imag_part() const;
    Value conjugate() const;
    Value magnitude() const;
    Value argument() const;
};

enum class TokenType {
    NUMBER, IDENTIFIER, PLUS, MINUS, STAR, SLASH, CARET, PERCENT,
    LPAREN, RPAREN, LBRACKET, RBRACKET, COMMA, SEMICOLON, BANG, EQUAL, LT, GT, LE, GE, NEQ,
    COLON_EQUAL,
    END_OF_INPUT, ERROR
};

struct Token {
    TokenType type;
    std::string text;
    int pos;

    Token() : type(TokenType::END_OF_INPUT), pos(0) {}
    Token(TokenType t, const std::string& s, int p) : type(t), text(s), pos(p) {}
};

class Lexer {
public:
    explicit Lexer(const std::string& input) : input_(input), pos_(0) {}

    std::vector<Token> tokenize();

private:
    std::string input_;
    int pos_;

    char peek() const;
    char advance();
    void skip_whitespace();
    Token read_number();
    Token read_identifier();
};

struct ASTNode;
using ASTPtr = std::shared_ptr<ASTNode>;

struct ASTNode {
    enum Type {
        NUMBER, CONSTANT, VARIABLE, BINOP, UNARYOP, FUNCTION, FACTORIAL, VEC_LITERAL, MAT_LITERAL
    };

    Type type;
    BigRat number;
    std::string name;
    char op;
    ASTPtr left;
    ASTPtr right;
    std::vector<ASTPtr> args;

    static ASTPtr make_num(const BigRat& n) {
        auto p = std::make_shared<ASTNode>();
        p->type = NUMBER;
        p->number = n;
        return p;
    }
    static ASTPtr make_const(const std::string& name) {
        auto p = std::make_shared<ASTNode>();
        p->type = CONSTANT;
        p->name = name;
        return p;
    }
    static ASTPtr make_var(const std::string& name) {
        auto p = std::make_shared<ASTNode>();
        p->type = VARIABLE;
        p->name = name;
        return p;
    }
    static ASTPtr make_binop(char op, ASTPtr left, ASTPtr right) {
        auto p = std::make_shared<ASTNode>();
        p->type = BINOP;
        p->op = op;
        p->left = std::move(left);
        p->right = std::move(right);
        return p;
    }
    static ASTPtr make_unaryop(char op, ASTPtr operand) {
        auto p = std::make_shared<ASTNode>();
        p->type = UNARYOP;
        p->op = op;
        p->left = std::move(operand);
        return p;
    }
    static ASTPtr make_func(const std::string& name, std::vector<ASTPtr> args) {
        auto p = std::make_shared<ASTNode>();
        p->type = FUNCTION;
        p->name = name;
        p->args = std::move(args);
        return p;
    }
    static ASTPtr make_factorial(ASTPtr operand) {
        auto p = std::make_shared<ASTNode>();
        p->type = FACTORIAL;
        p->left = std::move(operand);
        return p;
    }
    static ASTPtr make_vec_literal(std::vector<ASTPtr> components) {
        auto p = std::make_shared<ASTNode>();
        p->type = VEC_LITERAL;
        p->args = std::move(components);
        return p;
    }
    static ASTPtr make_mat_literal(std::vector<std::vector<ASTPtr>> rows) {
        auto p = std::make_shared<ASTNode>();
        p->type = MAT_LITERAL;
        p->mat_rows = std::move(rows);
        return p;
    }

    std::vector<std::vector<ASTPtr>> mat_rows;
};

class Parser {
public:
    explicit Parser(const std::vector<Token>& tokens) : tokens_(tokens), pos_(0) {}

    ASTPtr parse();

private:
    std::vector<Token> tokens_;
    size_t pos_;

    const Token& peek() const;
    Token advance();
    bool match(TokenType t);
    bool check(TokenType t) const;

    ASTPtr expression();
    ASTPtr comparison();
    ASTPtr additive();
    ASTPtr multiplicative();
    ASTPtr power();
    ASTPtr unary();
    ASTPtr postfix();
    ASTPtr primary();
    ASTPtr function_call(const std::string& name);
};

class Evaluator {
public:
    Evaluator() : angle_mode_(RAD), output_base_(10), last_ans_(BigRat(0)) {
        variables_["x"] = Value(BigRat(0));
        variables_["y"] = Value(BigRat(0));
        variables_["A"] = Value(BigRat(0));
        variables_["B"] = Value(BigRat(0));
        variables_["C"] = Value(BigRat(0));
        variables_["D"] = Value(BigRat(0));
        vec_variables_["VerA"] = Value::make_vector({});
        vec_variables_["VerB"] = Value::make_vector({});
        vec_variables_["VerC"] = Value::make_vector({});
        vec_variables_["VerD"] = Value::make_vector({});
    }

    enum AngleMode { DEG, RAD };

    Value evaluate(ASTPtr node);
    Value evaluate(const std::string& input);

    void set_angle_mode(AngleMode m) { angle_mode_ = m; }
    AngleMode angle_mode() const { return angle_mode_; }

    void set_variable(const std::string& name, const Value& v);
    Value get_variable(const std::string& name) const;

    void set_vec_variable(const std::string& name, const Value& v);
    Value get_vec_variable(const std::string& name) const;
    bool has_vec_variable(const std::string& name) const;

    void set_output_base(int b) { output_base_ = b; }
    int output_base() const { return output_base_; }

    Value get_last_ans() const { return last_ans_; }

    void set_user_function(const std::string& name, const std::vector<std::string>& params, ASTPtr body) {
        user_functions_[name] = {params, body};
    }
    bool has_user_function(const std::string& name) const {
        return user_functions_.find(name) != user_functions_.end();
    }
    std::pair<std::vector<std::string>, ASTPtr> get_user_function(const std::string& name) const {
        auto it = user_functions_.find(name);
        if (it != user_functions_.end()) return it->second;
        return {};
    }

    void set_latex_mode(bool m) { latex_mode_ = m; }
    bool latex_mode() const { return latex_mode_; }

    auto& user_functions() { return user_functions_; }

    std::string get_user_func_string(const std::string& name) const;

    static std::string format_result(const Value& v, int base = 10);
    static std::string format_latex(const Value& v);
    static std::string format_pretty_matrix(const Value& v);

private:
    AngleMode angle_mode_;
    int output_base_;
    bool latex_mode_ = false;
    std::map<std::string, Value> variables_;
    std::map<std::string, Value> vec_variables_;
    std::map<std::string, std::pair<std::vector<std::string>, ASTPtr>> user_functions_;
    Value last_ans_;

    Value eval_node(ASTPtr node);
    Value eval_binop(ASTPtr node);
    Value eval_unaryop(ASTPtr node);
    Value eval_function(ASTPtr node);
    Value eval_factorial(ASTPtr node);

    Value to_radians(const Value& v) const;
    Value from_radians(const Value& v) const;

    Value eval_sin(const Value& v);
    Value eval_cos(const Value& v);
    Value eval_tan(const Value& v);
    Value eval_asin(const Value& v);
    Value eval_acos(const Value& v);
    Value eval_atan(const Value& v);
    Value eval_sinh(const Value& v);
    Value eval_cosh(const Value& v);
    Value eval_tanh(const Value& v);
    Value eval_sqrt(const Value& v);
    Value eval_cbrt(const Value& v);
    Value eval_nrt(const Value& v, const Value& n);
    Value eval_abs(const Value& v);
    Value eval_ln(const Value& v);
    Value eval_lg(const Value& v);
    Value eval_log(const Value& base, const Value& v);
    Value eval_fact(const Value& v);
    Value eval_perm(const Value& n, const Value& k);
    Value eval_comb(const Value& n, const Value& k);
    Value eval_gcd(const Value& a, const Value& b);
    Value eval_lcm(const Value& a, const Value& b);
    Value eval_factor(const Value& v);
    Value eval_euler_phi(const Value& v);
    Value eval_int(ASTPtr node);
    Value eval_diff(ASTPtr node);
    Value eval_sum(ASTPtr node);
    Value eval_prod(ASTPtr node);
    Value eval_complex(const Value& r, const Value& i);
    Value eval_re(const Value& v);
    Value eval_im(const Value& v);
    Value eval_conj(const Value& v);
    Value eval_arg(const Value& v);
    Value eval_mod(const Value& v);
    Value eval_rand();
    Value eval_randint(const Value& a, const Value& b);
    Value eval_convert(const Value& v, const std::string& from, const std::string& to);

    Value eval_vecmod(const Value& v);
    Value eval_dot(const Value& a, const Value& b);
    Value eval_cross(const Value& a, const Value& b);
    Value eval_scalarmul(const Value& s, const Value& v);
    Value eval_mixed(const Value& a, const Value& b, const Value& c);
    Value eval_proj(const Value& a, const Value& b);
    Value eval_decompose(const Value& a, const Value& b, const Value& c);
    Value eval_decompose3d(const Value& a, const Value& b, const Value& c, const Value& d);

    Value eval_solve(const std::vector<ASTPtr>& args);
    Value eval_solve_quadratic(const Value& a, const Value& b, const Value& c);

    Value eval_matrix(const std::vector<ASTPtr>& args);
    Value eval_det(const Value& m);
    Value eval_inv(const Value& m);
    Value eval_eigen(const Value& m);
    Value eval_trace(const Value& m);
    Value eval_transpose(const Value& m);
    Value eval_identity(const Value& n);

    Value eval_mean(const std::vector<Value>& args);
    Value eval_stddev(const std::vector<Value>& args);
    Value eval_variance(const std::vector<Value>& args);
    Value eval_median(const std::vector<Value>& args);

    Value eval_simplify(ASTPtr node);

    ASTPtr sdiff(ASTPtr node, const std::string& var);
    Value eval_taylor(ASTPtr node);
    Value eval_limit(ASTPtr node);
    Value eval_inttable(ASTPtr node);
    Value eval_recur(ASTPtr node);
    Value eval_table(ASTPtr node);
    Value eval_lagrange(ASTPtr node);

    Value try_exact_trig(const BigRat& pi_coeff, int func);
    Value substitute(ASTPtr node, const std::string& var, const Value& val);
};

std::string run_repl();

#endif
