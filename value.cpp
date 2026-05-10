#include "ccalc.h"

ComplexVal::ComplexVal(const Value& r, const Value& i)
    : real(std::make_shared<Value>(r)), imag(std::make_shared<Value>(i)) {}

SurdsExpr::SurdsExpr() {
    terms.push_back(SurdsTerm(BigRat(0), BigInt(1)));
}

SurdsExpr::SurdsExpr(const BigRat& r) {
    if (!r.is_zero())
        terms.push_back(SurdsTerm(r, BigInt(1)));
}

SurdsExpr::SurdsExpr(const BigRat& c, const BigInt& radicand) {
    if (!c.is_zero()) {
        BigInt sf = make_square_free(radicand);
        terms.push_back(SurdsTerm(c, sf));
    }
}

BigInt SurdsExpr::make_square_free(const BigInt& n) {
    if (n.is_zero()) return BigInt(1);
    if (n < BigInt(0)) return n;
    if (n == BigInt(1)) return BigInt(1);
    BigInt result = n;
    for (int64_t i = 2; i * i <= 1000000; i++) {
        BigInt ii(i);
        BigInt ii2 = ii * ii;
        while (true) {
            auto dm = result.divmod(ii2);
            if (dm.second.is_zero()) {
                result = dm.first;
            } else {
                break;
            }
        }
    }
    return result;
}

void SurdsExpr::combine_like_terms() {
    std::map<std::string, BigRat> combined;
    for (auto& t : terms) {
        std::string key = t.radicand.to_string();
        if (combined.find(key) == combined.end())
            combined[key] = BigRat(0);
        combined[key] = combined[key] + t.coeff;
    }
    terms.clear();
    for (auto& [key, coeff] : combined) {
        if (!coeff.is_zero()) {
            terms.push_back(SurdsTerm(coeff, BigInt(key)));
        }
    }
    if (terms.empty()) {
        terms.push_back(SurdsTerm(BigRat(0), BigInt(1)));
    }
}

void SurdsExpr::simplify() {
    std::vector<SurdsTerm> new_terms;
    for (auto& t : terms) {
        if (t.coeff.is_zero()) continue;
        if (t.radicand == BigInt(1)) {
            new_terms.push_back(SurdsTerm(t.coeff, BigInt(1)));
            continue;
        }
        if (t.radicand.is_zero()) {
            new_terms.push_back(SurdsTerm(BigRat(0), BigInt(1)));
            continue;
        }
        if (t.radicand == BigInt(-1)) {
            new_terms.push_back(t);
            continue;
        }
        BigInt sf = make_square_free(t.radicand);
        if (sf == t.radicand) {
            new_terms.push_back(t);
        } else {
            BigInt quotient = t.radicand / sf;
            std::string qs = quotient.to_string();
            int64_t qv = std::stoll(qs);
            int64_t sq = (int64_t)std::sqrt((double)qv);
            if (sq * sq == qv) {
                new_terms.push_back(SurdsTerm(t.coeff * BigRat(sq), sf));
            } else {
                new_terms.push_back(SurdsTerm(t.coeff * BigRat(quotient), sf));
            }
        }
    }
    terms = new_terms;
    combine_like_terms();
}

SurdsExpr SurdsExpr::operator+(const SurdsExpr& o) const {
    SurdsExpr result;
    result.terms = terms;
    for (auto& t : o.terms) result.terms.push_back(t);
    result.combine_like_terms();
    return result;
}

SurdsExpr SurdsExpr::operator-(const SurdsExpr& o) const {
    SurdsExpr result;
    result.terms = terms;
    for (auto& t : o.terms) result.terms.push_back(SurdsTerm(-t.coeff, t.radicand));
    result.combine_like_terms();
    return result;
}

SurdsExpr SurdsExpr::operator*(const SurdsExpr& o) const {
    bool has_const_a = false, has_const_b = false;
    bool has_e_a = false, has_e_b = false;
    for (auto& t : terms) {
        if (t.radicand == BigInt(-1)) has_const_a = true;
        if (t.radicand == BigInt(-2)) has_e_a = true;
    }
    for (auto& t : o.terms) {
        if (t.radicand == BigInt(-1)) has_const_b = true;
        if (t.radicand == BigInt(-2)) has_e_b = true;
    }
    if ((has_const_a && has_const_b) || (has_e_a && has_e_b) ||
        (has_const_a && has_e_b) || (has_e_a && has_const_b)) {
        BigFloat fa = to_float();
        BigFloat fb = o.to_float();
        BigFloat r = fa * fb;
        return SurdsExpr(BigRat(0));
    }
    SurdsExpr result;
    result.terms.clear();
    for (auto& a : terms) {
        for (auto& b : o.terms) {
            BigRat coeff = a.coeff * b.coeff;
            bool a_is_const = (a.radicand == BigInt(-1) || a.radicand == BigInt(-2));
            bool b_is_const = (b.radicand == BigInt(-1) || b.radicand == BigInt(-2));
            if (a_is_const || b_is_const) {
                if (a_is_const && b_is_const) continue;
                BigInt const_radicand = a_is_const ? a.radicand : b.radicand;
                BigInt other_radicand = a_is_const ? b.radicand : a.radicand;
                if (other_radicand == BigInt(1)) {
                    result.terms.push_back(SurdsTerm(coeff, const_radicand));
                } else {
                    BigFloat fa = to_float();
                    BigFloat fb = o.to_float();
                    BigFloat r = fa * fb;
                    return SurdsExpr(BigRat(0));
                }
            } else {
                BigInt radicand = a.radicand * b.radicand;
                if (!coeff.is_zero()) {
                    result.terms.push_back(SurdsTerm(coeff, radicand));
                }
            }
        }
    }
    result.simplify();
    return result;
}

SurdsExpr SurdsExpr::operator/(const SurdsExpr& o) const {
    if (o.is_zero()) throw CalcError("Division by zero");
    if (o.terms.size() == 1 && o.terms[0].radicand == BigInt(1)) {
        BigRat divisor = o.terms[0].coeff;
        SurdsExpr result;
        result.terms.clear();
        for (auto& t : terms) {
            result.terms.push_back(SurdsTerm(t.coeff / divisor, t.radicand));
        }
        return result;
    }
    BigFloat a = to_float();
    BigFloat b = o.to_float();
    return SurdsExpr(BigRat(1)) * SurdsExpr(BigRat(0), BigInt(1));
}

SurdsExpr SurdsExpr::operator-() const {
    SurdsExpr result;
    result.terms.clear();
    for (auto& t : terms) {
        result.terms.push_back(SurdsTerm(-t.coeff, t.radicand));
    }
    return result;
}

bool SurdsExpr::operator==(const SurdsExpr& o) const {
    SurdsExpr diff = *this - o;
    return diff.is_zero();
}

bool SurdsExpr::is_zero() const {
    for (auto& t : terms) {
        if (!t.coeff.is_zero()) return false;
    }
    return true;
}

bool SurdsExpr::is_rational() const {
    for (auto& t : terms) {
        if (t.radicand != BigInt(1)) return false;
    }
    return true;
}

bool SurdsExpr::is_positive() const {
    return to_float() > BigFloat(0);
}

bool SurdsExpr::is_negative() const {
    return to_float() < BigFloat(0);
}

BigRat SurdsExpr::to_rational() const {
    BigRat result(0);
    for (auto& t : terms) {
        if (t.radicand == BigInt(1)) {
            result = result + t.coeff;
        } else {
            throw CalcError("Cannot convert surds expression to rational");
        }
    }
    return result;
}

BigFloat SurdsExpr::to_float(int prec) const {
    int p = prec > 0 ? prec : g_precision;
    BigFloat result(0);
    for (auto& t : terms) {
        BigFloat c(t.coeff);
        if (t.radicand == BigInt(1)) {
            result = result + c;
        } else if (t.radicand == BigInt(-1)) {
            result = result + c * BigFloat::pi(p);
        } else if (t.radicand == BigInt(-2)) {
            result = result + c * BigFloat::e_val(p);
        } else {
            BigFloat r(t.radicand);
            BigFloat sr = BigFloat::sqrt_val(r, p);
            result = result + c * sr;
        }
    }
    return result;
}

std::string SurdsExpr::to_string() const {
    if (is_zero()) return "0";
    std::string result;
    bool first = true;
    auto sorted_terms = terms;
    std::sort(sorted_terms.begin(), sorted_terms.end(),
              [](const SurdsTerm& a, const SurdsTerm& b) {
                  return a.radicand > b.radicand;
              });
    for (auto& t : sorted_terms) {
        if (t.coeff.is_zero()) continue;
        bool is_rat = t.radicand == BigInt(1);
        bool is_pi = t.radicand == BigInt(-1);
        bool is_e = t.radicand == BigInt(-2);
        BigRat abs_coeff = t.coeff.abs();
        bool neg = t.coeff < BigRat(0);
        if (first) {
            if (neg) result += "-";
            first = false;
        } else {
            result += neg ? " - " : " + ";
        }
        if (is_rat) {
            result += abs_coeff.to_string();
        } else {
            BigInt p = abs_coeff.numerator().abs();
            BigInt q = abs_coeff.denominator();
            bool p_is_one = p == BigInt(1);
            bool q_is_one = q == BigInt(1);
            std::string sym;
            if (is_pi) sym = "pi";
            else if (is_e) sym = "e";
            else sym = "sqrt(" + t.radicand.to_string() + ")";
            if (p_is_one && q_is_one) {
                result += sym;
            } else if (p_is_one && !q_is_one) {
                result += sym + "/" + q.to_string();
            } else if (!p_is_one && q_is_one) {
                result += p.to_string() + "*" + sym;
            } else {
                result += p.to_string() + "*" + sym + "/" + q.to_string();
            }
        }
    }
    if (result.empty()) return "0";
    return result;
}

Value Value::make_complex(const Value& r, const Value& i) {
    Value v;
    v.type = COMPLEX;
    v.complex = ComplexVal(r, i);
    return v;
}

bool Value::is_exact() const {
    if (type == SURDS) return true;
    if (type == COMPLEX) {
        return complex.real->is_exact() && complex.imag->is_exact();
    }
    return false;
}

bool Value::is_rational() const {
    if (type == SURDS) return surds.is_rational();
    return false;
}

bool Value::is_zero() const {
    if (type == SURDS) return surds.is_zero();
    if (type == FLOAT) return float_val.is_zero();
    if (type == COMPLEX)
        return complex.real->is_zero() && complex.imag->is_zero();
    if (type == VECTOR) {
        for (auto& c : vec) if (!c.is_zero()) return false;
        return true;
    }
    if (type == MATRIX) {
        for (auto& row : mat) for (auto& c : row) if (!c.is_zero()) return false;
        return true;
    }
    return false;
}

BigRat Value::to_rational() const {
    if (type == SURDS) return surds.to_rational();
    throw CalcError("Cannot convert to rational");
}

BigFloat Value::to_float(int prec) const {
    int p = prec > 0 ? prec : g_precision;
    if (type == SURDS) return surds.to_float(p);
    if (type == FLOAT) return float_val;
    if (type == COMPLEX)
        throw CalcError("Cannot convert complex to float");
    throw CalcError("Cannot convert error to float");
}

std::string Value::to_string() const {
    if (type == ERROR) return "Error: " + error_msg;
    if (type == STRING) return error_msg;
    if (type == VECTOR) {
        std::string r = "(";
        for (size_t i = 0; i < vec.size(); i++) {
            if (i > 0) r += ", ";
            r += vec[i].to_string();
        }
        r += ")";
        return r;
    }
    if (type == MATRIX) {
        std::string r = "[";
        for (size_t i = 0; i < mat.size(); i++) {
            if (i > 0) r += "; ";
            r += "[";
            for (size_t j = 0; j < mat[i].size(); j++) {
                if (j > 0) r += ", ";
                r += mat[i][j].to_string();
            }
            r += "]";
        }
        r += "]";
        return r;
    }
    if (type == SURDS) {
        if (surds.is_rational()) {
            BigRat r = surds.to_rational();
            if (r.denominator() == BigInt(1)) return r.numerator().to_string();
            return r.to_string();
        }
        return surds.to_string();
    }
    if (type == FLOAT) return float_val.to_string();
    if (type == COMPLEX) {
        std::string r = complex.real->to_string();
        std::string i = complex.imag->to_string();
        if (complex.imag->is_zero()) return r;
        if (complex.real->is_zero()) {
            if (i == "1") return "i";
            if (i == "-1") return "-i";
            return i + "*i";
        }
        if (i == "1") return r + " + i";
        if (i == "-1") return r + " - i";
        if (complex.imag->to_string()[0] == '-') {
            return r + " - " + i.substr(1) + "*i";
        }
        return r + " + " + i + "*i";
    }
    return "undefined";
}

Value Value::operator+(const Value& o) const {
    if (is_error()) return *this;
    if (o.is_error()) return o;
    if (type == VECTOR || o.type == VECTOR) {
        if (type == VECTOR && o.type == VECTOR) {
            if (vec_dim() != o.vec_dim())
                return make_error("Vector dimension mismatch in addition");
            std::vector<Value> result;
            for (int i = 0; i < vec_dim(); i++)
                result.push_back(vec[i] + o.vec[i]);
            return make_vector(result);
        }
        return make_error("Cannot add vector and scalar");
    }
    if (type == MATRIX || o.type == MATRIX) {
        if (type == MATRIX && o.type == MATRIX) {
            if (mat_rows() != o.mat_rows() || mat_cols() != o.mat_cols())
                return make_error("Matrix dimension mismatch in addition");
            std::vector<std::vector<Value>> result;
            for (int i = 0; i < mat_rows(); i++) {
                std::vector<Value> row;
                for (int j = 0; j < mat_cols(); j++)
                    row.push_back(mat[i][j] + o.mat[i][j]);
                result.push_back(row);
            }
            return make_matrix(result);
        }
        return make_error("Cannot add matrix and scalar");
    }
    if (type == COMPLEX || o.type == COMPLEX) {
        Value r1 = type == COMPLEX ? *complex.real : *this;
        Value i1 = type == COMPLEX ? *complex.imag : Value(BigRat(0));
        Value r2 = o.type == COMPLEX ? *o.complex.real : o;
        Value i2 = o.type == COMPLEX ? *o.complex.imag : Value(BigRat(0));
        return make_complex(r1 + r2, i1 + i2);
    }
    if (type == SURDS && o.type == SURDS) {
        SurdsExpr result = surds + o.surds;
        if (result.is_rational()) return Value(result.to_rational());
        return Value(result);
    }
    if (type == FLOAT || o.type == FLOAT) {
        return Value(to_float() + o.to_float());
    }
    return Value(BigRat(0));
}

Value Value::operator-(const Value& o) const {
    if (is_error()) return *this;
    if (o.is_error()) return o;
    if (type == VECTOR || o.type == VECTOR) {
        if (type == VECTOR && o.type == VECTOR) {
            if (vec_dim() != o.vec_dim())
                return make_error("Vector dimension mismatch in subtraction");
            std::vector<Value> result;
            for (int i = 0; i < vec_dim(); i++)
                result.push_back(vec[i] - o.vec[i]);
            return make_vector(result);
        }
        return make_error("Cannot subtract vector and scalar");
    }
    if (type == MATRIX || o.type == MATRIX) {
        if (type == MATRIX && o.type == MATRIX) {
            if (mat_rows() != o.mat_rows() || mat_cols() != o.mat_cols())
                return make_error("Matrix dimension mismatch in subtraction");
            std::vector<std::vector<Value>> result;
            for (int i = 0; i < mat_rows(); i++) {
                std::vector<Value> row;
                for (int j = 0; j < mat_cols(); j++)
                    row.push_back(mat[i][j] - o.mat[i][j]);
                result.push_back(row);
            }
            return make_matrix(result);
        }
        return make_error("Cannot subtract matrix and scalar");
    }
    if (type == COMPLEX || o.type == COMPLEX) {
        Value r1 = type == COMPLEX ? *complex.real : *this;
        Value i1 = type == COMPLEX ? *complex.imag : Value(BigRat(0));
        Value r2 = o.type == COMPLEX ? *o.complex.real : o;
        Value i2 = o.type == COMPLEX ? *o.complex.imag : Value(BigRat(0));
        return make_complex(r1 - r2, i1 - i2);
    }
    if (type == SURDS && o.type == SURDS) {
        SurdsExpr result = surds - o.surds;
        if (result.is_rational()) return Value(result.to_rational());
        return Value(result);
    }
    if (type == FLOAT || o.type == FLOAT) {
        return Value(to_float() - o.to_float());
    }
    return Value(BigRat(0));
}

Value Value::operator*(const Value& o) const {
    if (is_error()) return *this;
    if (o.is_error()) return o;
    if (type == VECTOR || o.type == VECTOR) {
        if (type == VECTOR && o.type == VECTOR)
            return make_error("Use dot(a,b) for dot product");
        if (type == VECTOR) {
            std::vector<Value> result;
            for (auto& c : vec) result.push_back(c * o);
            return make_vector(result);
        }
        std::vector<Value> result;
        for (auto& c : o.vec) result.push_back(*this * c);
        return make_vector(result);
    }
    if (type == MATRIX && o.type == MATRIX) {
        if (mat_cols() != o.mat_rows())
            return make_error("Matrix dimension mismatch in multiplication");
        std::vector<std::vector<Value>> result;
        for (int i = 0; i < mat_rows(); i++) {
            std::vector<Value> row;
            for (int j = 0; j < o.mat_cols(); j++) {
                Value sum(BigRat(0));
                for (int k = 0; k < mat_cols(); k++)
                    sum = sum + mat[i][k] * o.mat[k][j];
                row.push_back(sum);
            }
            result.push_back(row);
        }
        return make_matrix(result);
    }
    if (type == MATRIX) {
        std::vector<std::vector<Value>> result;
        for (auto& row : mat) {
            std::vector<Value> nr;
            for (auto& c : row) nr.push_back(c * o);
            result.push_back(nr);
        }
        return make_matrix(result);
    }
    if (o.type == MATRIX) {
        std::vector<std::vector<Value>> result;
        for (auto& row : o.mat) {
            std::vector<Value> nr;
            for (auto& c : row) nr.push_back(*this * c);
            result.push_back(nr);
        }
        return make_matrix(result);
    }
    if (type == COMPLEX || o.type == COMPLEX) {
        Value r1 = type == COMPLEX ? *complex.real : *this;
        Value i1 = type == COMPLEX ? *complex.imag : Value(BigRat(0));
        Value r2 = o.type == COMPLEX ? *o.complex.real : o;
        Value i2 = o.type == COMPLEX ? *o.complex.imag : Value(BigRat(0));
        Value nr = r1 * r2 - i1 * i2;
        Value ni = r1 * i2 + i1 * r2;
        return make_complex(nr, ni);
    }
    if (type == SURDS && o.type == SURDS) {
        SurdsExpr result = surds * o.surds;
        if (result.is_rational()) return Value(result.to_rational());
        return Value(result);
    }
    if (type == FLOAT || o.type == FLOAT) {
        return Value(to_float() * o.to_float());
    }
    return Value(BigRat(0));
}

Value Value::operator/(const Value& o) const {
    if (is_error()) return *this;
    if (o.is_error()) return o;
    if (o.is_zero()) return make_error("Division by zero");
    if (type == COMPLEX || o.type == COMPLEX) {
        Value r1 = type == COMPLEX ? *complex.real : *this;
        Value i1 = type == COMPLEX ? *complex.imag : Value(BigRat(0));
        Value r2 = o.type == COMPLEX ? *o.complex.real : o;
        Value i2 = o.type == COMPLEX ? *o.complex.imag : Value(BigRat(0));
        Value denom = r2 * r2 + i2 * i2;
        Value nr = (r1 * r2 + i1 * i2) / denom;
        Value ni = (i1 * r2 - r1 * i2) / denom;
        return make_complex(nr, ni);
    }
    if (type == SURDS && o.type == SURDS && o.surds.is_rational()) {
        BigRat divisor = o.surds.to_rational();
        SurdsExpr result;
        result.terms.clear();
        for (auto& t : surds.terms) {
            result.terms.push_back(SurdsTerm(t.coeff / divisor, t.radicand));
        }
        result.simplify();
        if (result.is_rational()) return Value(result.to_rational());
        return Value(result);
    }
    return Value(to_float() / o.to_float());
}

Value Value::operator-() const {
    if (is_error()) return *this;
    if (type == VECTOR) {
        std::vector<Value> result;
        for (auto& c : vec) result.push_back(-c);
        return make_vector(result);
    }
    if (type == COMPLEX) return make_complex(-*complex.real, -*complex.imag);
    if (type == SURDS) return Value(-surds);
    if (type == FLOAT) return Value(-float_val);
    return *this;
}

Value Value::real_part() const {
    if (type == COMPLEX) return *complex.real;
    return *this;
}

Value Value::imag_part() const {
    if (type == COMPLEX) return *complex.imag;
    return Value(BigRat(0));
}

Value Value::conjugate() const {
    if (type == COMPLEX) return make_complex(*complex.real, -*complex.imag);
    return *this;
}

Value Value::magnitude() const {
    if (type == COMPLEX) {
        Value r = *complex.real;
        Value i = *complex.imag;
        return Value(BigFloat::sqrt_val(
            (r.to_float() * r.to_float() + i.to_float() * i.to_float())));
    }
    if (type == SURDS) {
        if (surds.is_rational()) {
            BigRat r = surds.to_rational();
            if (r < BigRat(0)) return Value(-r);
            return Value(r);
        }
        return Value(surds.to_float().abs());
    }
    return Value(float_val.abs());
}

Value Value::argument() const {
    if (type == COMPLEX) {
        return Value(BigFloat::atan_val(
            complex.imag->to_float() / complex.real->to_float()));
    }
    if (type == SURDS) {
        if (surds.is_rational()) {
            BigRat r = surds.to_rational();
            if (r < BigRat(0)) return Value(SurdsExpr(BigRat(-1), BigInt(-1)));
            return Value(BigRat(0));
        }
        if (surds.is_negative()) return Value(SurdsExpr(BigRat(1), BigInt(-1)));
        return Value(BigRat(0));
    }
    if (float_val < BigFloat(0)) return Value(BigFloat::pi());
    return Value(BigRat(0));
}
