#include "ccalc.h"

int g_precision = 110;

BigInt::BigInt(int64_t n) : neg_(n < 0) {
    if (n < 0) n = -n;
    if (n == 0) {
        digits_.push_back(0);
    } else {
        while (n > 0) {
            digits_.push_back(n % BASE);
            n /= BASE;
        }
    }
}

BigInt::BigInt(const std::string& s) : neg_(false) {
    int start = 0;
    if (start < (int)s.size() && s[start] == '-') {
        neg_ = true;
        start++;
    } else if (start < (int)s.size() && s[start] == '+') {
        start++;
    }
    while (start < (int)s.size() && s[start] == '0') start++;
    if (start >= (int)s.size() || (start == (int)s.size() - 1 && s[start] == '.')) {
        digits_.push_back(0);
        neg_ = false;
        return;
    }
    std::string digits_str = s.substr(start);
    std::string clean;
    for (char c : digits_str) {
        if (c == '.') break;
        if (c >= '0' && c <= '9') clean += c;
    }
    if (clean.empty()) {
        digits_.push_back(0);
        neg_ = false;
        return;
    }
    for (int i = (int)clean.size(); i > 0; i -= BASE_DIGITS) {
        int begin = std::max(0, i - BASE_DIGITS);
        digits_.push_back(std::stoll(clean.substr(begin, i - begin)));
    }
    trim();
}

void BigInt::trim() {
    while (digits_.size() > 1 && digits_.back() == 0)
        digits_.pop_back();
    if (digits_.size() == 1 && digits_[0] == 0)
        neg_ = false;
}

BigInt BigInt::abs() const {
    BigInt r = *this;
    r.neg_ = false;
    return r;
}

int BigInt::sign() const {
    if (is_zero()) return 0;
    return neg_ ? -1 : 1;
}

bool BigInt::is_zero() const {
    return digits_.size() == 1 && digits_[0] == 0;
}

BigInt BigInt::operator-() const {
    BigInt r = *this;
    if (!r.is_zero()) r.neg_ = !r.neg_;
    return r;
}

int BigInt::raw_compare(const BigInt& o) const {
    if (digits_.size() != o.digits_.size())
        return (int)digits_.size() - (int)o.digits_.size();
    for (int i = (int)digits_.size() - 1; i >= 0; i--) {
        if (digits_[i] != o.digits_[i])
            return digits_[i] > o.digits_[i] ? 1 : -1;
    }
    return 0;
}

BigInt BigInt::raw_add(const BigInt& o) const {
    BigInt res;
    res.neg_ = false;
    res.digits_.clear();
    int64_t carry = 0;
    size_t n = std::max(digits_.size(), o.digits_.size());
    for (size_t i = 0; i < n || carry; i++) {
        int64_t sum = carry;
        if (i < digits_.size()) sum += digits_[i];
        if (i < o.digits_.size()) sum += o.digits_[i];
        res.digits_.push_back(sum % BASE);
        carry = sum / BASE;
    }
    res.trim();
    return res;
}

BigInt BigInt::raw_sub(const BigInt& o) const {
    BigInt res;
    res.neg_ = false;
    res.digits_.clear();
    int64_t borrow = 0;
    for (size_t i = 0; i < digits_.size(); i++) {
        int64_t diff = digits_[i] - borrow;
        if (i < o.digits_.size()) diff -= o.digits_[i];
        if (diff < 0) {
            diff += BASE;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res.digits_.push_back(diff);
    }
    res.trim();
    return res;
}

BigInt BigInt::operator+(const BigInt& o) const {
    if (neg_ == o.neg_) {
        BigInt r = raw_add(o);
        r.neg_ = neg_;
        return r;
    }
    int cmp = raw_compare(o);
    if (cmp == 0) return BigInt(0);
    if (cmp > 0) {
        BigInt r = raw_sub(o);
        r.neg_ = neg_;
        return r;
    } else {
        BigInt r = o.raw_sub(*this);
        r.neg_ = o.neg_;
        return r;
    }
}

BigInt BigInt::operator-(const BigInt& o) const {
    return *this + (-o);
}

BigInt BigInt::operator*(const BigInt& o) const {
    if (is_zero() || o.is_zero()) return BigInt(0);
    BigInt res;
    res.neg_ = neg_ != o.neg_;
    const BigInt& a = *this;
    const BigInt& b = o;
    int n = (int)std::max(a.digits_.size(), b.digits_.size());
    if (n < 64) {
        res.digits_.assign(a.digits_.size() + b.digits_.size(), 0);
        for (size_t i = 0; i < a.digits_.size(); i++) {
            int64_t carry = 0;
            for (size_t j = 0; j < b.digits_.size() || carry; j++) {
                int64_t cur = res.digits_[i + j] + carry;
                if (j < b.digits_.size()) cur += a.digits_[i] * b.digits_[j];
                res.digits_[i + j] = cur % BASE;
                carry = cur / BASE;
            }
        }
        res.trim();
        return res;
    }
    res.digits_ = karatsuba(a.digits_, b.digits_);
    res.trim();
    return res;
}

std::vector<int64_t> BigInt::karatsuba(const std::vector<int64_t>& x,
                                        const std::vector<int64_t>& y) {
    size_t n = std::max(x.size(), y.size());
    if (n < 32) {
        std::vector<int64_t> res(x.size() + y.size(), 0);
        for (size_t i = 0; i < x.size(); i++) {
            int64_t carry = 0;
            for (size_t j = 0; j < y.size() || carry; j++) {
                int64_t cur = res[i + j] + carry;
                if (j < y.size()) cur += x[i] * y[j];
                res[i + j] = cur % BASE;
                carry = cur / BASE;
            }
        }
        while (res.size() > 1 && res.back() == 0) res.pop_back();
        return res;
    }
    size_t m = n / 2;
    std::vector<int64_t> x0(x.begin(), x.begin() + std::min(m, x.size()));
    std::vector<int64_t> y0(y.begin(), y.begin() + std::min(m, y.size()));
    std::vector<int64_t> x1(x.begin() + std::min(m, x.size()), x.end());
    std::vector<int64_t> y1(y.begin() + std::min(m, y.size()), y.end());
    if (x0.empty()) x0.push_back(0);
    if (y0.empty()) y0.push_back(0);
    if (x1.empty()) x1.push_back(0);
    if (y1.empty()) y1.push_back(0);
    auto z0 = karatsuba(x0, y0);
    auto z2 = karatsuba(x1, y1);
    auto xs = add_vec(x0, x1);
    auto ys = add_vec(y0, y1);
    auto z1_full = karatsuba(xs, ys);
    auto z1 = sub_vec(sub_vec(z1_full, z0), z2);
    auto result = add_vec(z0, shift(z1, m));
    result = add_vec(result, shift(z2, 2 * m));
    return result;
}

std::vector<int64_t> BigInt::add_vec(const std::vector<int64_t>& a,
                                      const std::vector<int64_t>& b) {
    std::vector<int64_t> res;
    int64_t carry = 0;
    size_t i = 0;
    while (i < a.size() || i < b.size() || carry) {
        int64_t sum = carry;
        if (i < a.size()) sum += a[i];
        if (i < b.size()) sum += b[i];
        res.push_back(sum % BASE);
        carry = sum / BASE;
        i++;
    }
    return res;
}

std::vector<int64_t> BigInt::sub_vec(const std::vector<int64_t>& a,
                                      const std::vector<int64_t>& b) {
    std::vector<int64_t> res;
    int64_t borrow = 0;
    for (size_t i = 0; i < a.size(); i++) {
        int64_t diff = a[i] - borrow;
        if (i < b.size()) diff -= b[i];
        if (diff < 0) {
            diff += BASE;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res.push_back(diff);
    }
    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

std::vector<int64_t> BigInt::shift(const std::vector<int64_t>& v, size_t offset) {
    std::vector<int64_t> res(offset, 0);
    res.insert(res.end(), v.begin(), v.end());
    return res;
}

std::pair<BigInt, BigInt> BigInt::divmod(const BigInt& o) const {
    if (o.is_zero()) throw CalcError("Division by zero");
    BigInt a = abs();
    BigInt b = o.abs();
    if (a < b) {
        BigInt q(0);
        BigInt r = *this;
        return {q, r};
    }
    BigInt q;
    q.digits_.assign(a.digits_.size(), 0);
    BigInt r(0);
    for (int i = (int)a.digits_.size() - 1; i >= 0; i--) {
        r.digits_.insert(r.digits_.begin(), a.digits_[i]);
        r.trim();
        int64_t lo = 0, hi = BASE - 1;
        while (lo <= hi) {
            int64_t mid = (lo + hi) / 2;
            BigInt t = b * BigInt(mid);
            if (t <= r) {
                q.digits_[i] = mid;
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        r = r - b * BigInt(q.digits_[i]);
    }
    q.neg_ = neg_ != o.neg_;
    r.neg_ = neg_;
    q.trim();
    r.trim();
    if (r.neg_ && !r.is_zero()) {
        r = r + o.abs();
        q = q - BigInt(1);
    }
    return {q, r};
}

BigInt BigInt::operator/(const BigInt& o) const {
    return divmod(o).first;
}

BigInt BigInt::operator%(const BigInt& o) const {
    return divmod(o).second;
}

bool BigInt::operator==(const BigInt& o) const {
    return neg_ == o.neg_ && raw_compare(o) == 0;
}

bool BigInt::operator<(const BigInt& o) const {
    if (neg_ != o.neg_) return neg_;
    int cmp = raw_compare(o);
    return neg_ ? cmp > 0 : cmp < 0;
}

std::string BigInt::to_string() const {
    if (is_zero()) return "0";
    std::string s;
    if (neg_) s = "-";
    s += std::to_string(digits_.back());
    for (int i = (int)digits_.size() - 2; i >= 0; i--) {
        std::string t = std::to_string(digits_[i]);
        s += std::string(BASE_DIGITS - t.size(), '0') + t;
    }
    return s;
}

std::string BigInt::to_decimal_string() const {
    return to_string();
}

std::string BigInt::to_base_string(int base) const {
    if (base < 2 || base > 36) return to_string();
    if (is_zero()) return "0";
    const char digits[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    BigInt n = abs();
    std::string result;
    while (!n.is_zero()) {
        auto dm = n.divmod(BigInt(base));
        int d = dm.second.to_int64();
        result = digits[d] + result;
        n = dm.first;
    }
    if (neg_) result = "-" + result;
    return result;
}

BigInt BigInt::from_base_string(const std::string& s, int base) {
    if (base < 2 || base > 36) return BigInt(s);
    if (s.empty()) return BigInt(0);
    bool neg = false;
    size_t start = 0;
    if (s[0] == '-') { neg = true; start = 1; }
    else if (s[0] == '+') { start = 1; }
    BigInt result(0);
    for (size_t i = start; i < s.size(); i++) {
        char c = std::toupper(s[i]);
        int d = 0;
        if (c >= '0' && c <= '9') d = c - '0';
        else if (c >= 'A' && c <= 'Z') d = c - 'A' + 10;
        else continue;
        if (d >= base) continue;
        result = result * BigInt(base) + BigInt(d);
    }
    if (neg) result = -result;
    return result;
}

int64_t BigInt::to_int64() const {
    if (is_zero()) return 0;
    int64_t r = 0;
    for (int i = (int)digits_.size() - 1; i >= 0; i--) {
        r = r * BASE + digits_[i];
    }
    return neg_ ? -r : r;
}

BigInt BigInt::gcd(BigInt a, BigInt b) {
    a = a.abs();
    b = b.abs();
    while (!b.is_zero()) {
        BigInt t = b;
        b = a % b;
        a = t;
    }
    return a;
}

BigInt BigInt::lcm(const BigInt& a, const BigInt& b) {
    if (a.is_zero() || b.is_zero()) return BigInt(0);
    return a.abs() / gcd(a, b) * b.abs();
}

BigInt BigInt::pow(BigInt base, BigInt exp) {
    if (exp.is_zero()) return BigInt(1);
    if (exp < BigInt(0)) throw CalcError("Negative exponent in integer power");
    BigInt result(1);
    while (!exp.is_zero()) {
        if (exp.digits_[0] % 2 == 1) {
            result = result * base;
        }
        base = base * base;
        exp = exp / BigInt(2);
    }
    return result;
}

BigInt BigInt::factorial(int64_t n) {
    if (n < 0) throw CalcError("Factorial of negative number");
    BigInt result(1);
    for (int64_t i = 2; i <= n; i++) {
        result = result * BigInt(i);
    }
    return result;
}

BigRat::BigRat(int64_t n, int64_t d) : num_(n), den_(d) {
    if (den_.is_zero()) throw CalcError("Division by zero: denominator is 0");
    if (den_ < BigInt(0)) {
        num_ = -num_;
        den_ = -den_;
    }
    simplify();
}

BigRat::BigRat(const BigInt& n, const BigInt& d) : num_(n), den_(d) {
    if (den_.is_zero()) throw CalcError("Division by zero: denominator is 0");
    if (den_ < BigInt(0)) {
        num_ = -num_;
        den_ = -den_;
    }
    simplify();
}

BigRat::BigRat(const std::string& s) : num_(0), den_(1) {
    size_t slash = s.find('/');
    if (slash == std::string::npos) {
        size_t dot = s.find('.');
        if (dot == std::string::npos) {
            num_ = BigInt(s);
            den_ = BigInt(1);
        } else {
            std::string int_part = s.substr(0, dot);
            std::string frac_part = s.substr(dot + 1);
            BigInt numerator(int_part + frac_part);
            BigInt denominator(1);
            for (size_t i = 0; i < frac_part.size(); i++)
                denominator = denominator * BigInt(10);
            num_ = numerator;
            den_ = denominator;
            simplify();
        }
    } else {
        BigRat left(s.substr(0, slash));
        BigRat right(s.substr(slash + 1));
        *this = left / right;
    }
}

void BigRat::simplify() {
    if (num_.is_zero()) {
        den_ = BigInt(1);
        return;
    }
    BigInt g = BigInt::gcd(num_.abs(), den_.abs());
    num_ = num_ / g;
    den_ = den_ / g;
    if (den_ < BigInt(0)) {
        num_ = -num_;
        den_ = -den_;
    }
}

BigRat BigRat::operator+(const BigRat& o) const {
    return BigRat(num_ * o.den_ + o.num_ * den_, den_ * o.den_);
}

BigRat BigRat::operator-(const BigRat& o) const {
    return BigRat(num_ * o.den_ - o.num_ * den_, den_ * o.den_);
}

BigRat BigRat::operator*(const BigRat& o) const {
    return BigRat(num_ * o.num_, den_ * o.den_);
}

BigRat BigRat::operator/(const BigRat& o) const {
    if (o.num_.is_zero()) throw CalcError("Division by zero");
    return BigRat(num_ * o.den_, den_ * o.num_);
}

BigRat BigRat::operator-() const {
    BigRat r;
    r.num_ = -num_;
    r.den_ = den_;
    return r;
}

bool BigRat::operator==(const BigRat& o) const {
    return num_ == o.num_ && den_ == o.den_;
}

bool BigRat::operator<(const BigRat& o) const {
    return num_ * o.den_ < o.num_ * den_;
}

BigRat BigRat::abs() const {
    BigRat r;
    r.num_ = num_.abs();
    r.den_ = den_;
    return r;
}

int BigRat::sign() const {
    return num_.sign();
}

bool BigRat::is_zero() const {
    return num_.is_zero();
}

std::string BigRat::to_string() const {
    if (den_ == BigInt(1)) return num_.to_string();
    return num_.to_string() + "/" + den_.to_string();
}

std::string BigRat::to_mixed_string() const {
    if (den_ == BigInt(1)) return num_.to_string();
    BigInt q = num_ / den_;
    BigInt r = num_ % den_;
    if (r.is_zero()) return q.to_string();
    if (q.is_zero()) return r.abs().to_string() + "/" + den_.to_string();
    std::string s;
    if (num_ < BigInt(0)) s += "-";
    s += q.abs().to_string() + "(" + r.abs().to_string() + "/" + den_.to_string() + ")";
    return s;
}

std::string BigRat::to_decimal_string(int precision) const {
    int p = precision > 0 ? precision : g_precision;
    if (is_zero()) return "0";
    BigInt n = num_.abs();
    BigInt d = den_;
    BigInt int_part = n / d;
    BigInt rem = n % d;
    std::string s;
    if (num_ < BigInt(0)) s += "-";
    s += int_part.to_string();
    if (rem.is_zero()) return s;
    s += ".";
    for (int i = 0; i < p + 5; i++) {
        rem = rem * BigInt(10);
        BigInt digit = rem / d;
        rem = rem % d;
        s += digit.to_string();
        if (rem.is_zero()) break;
    }
    if ((int)s.size() - (int)s.find('.') - 1 > p) {
        s = s.substr(0, s.find('.') + p + 1);
    }
    return s;
}

BigRat BigRat::pow(const BigRat& base, int64_t exp) {
    if (exp == 0) return BigRat(1);
    if (exp < 0) return BigRat(1) / pow(base, -exp);
    BigRat result(1);
    BigRat b = base;
    while (exp > 0) {
        if (exp % 2 == 1) result = result * b;
        b = b * b;
        exp /= 2;
    }
    return result;
}

BigFloat::BigFloat(int64_t n) : mantissa_(n), exp_(0) {
    normalize();
}

BigFloat::BigFloat(double n) {
    if (n == 0.0) {
        mantissa_ = BigInt(0);
        exp_ = 0;
        return;
    }
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(17) << n;
    std::string s = oss.str();
    size_t e_pos = s.find('e');
    if (e_pos == std::string::npos) e_pos = s.find('E');
    std::string mantissa_str = s.substr(0, e_pos);
    int64_t e = std::stoll(s.substr(e_pos + 1));
    mantissa_str.erase(std::remove(mantissa_str.begin(), mantissa_str.end(), '.'), mantissa_str.end());
    mantissa_ = BigInt(mantissa_str);
    exp_ = e - 17 + 1;
    normalize();
}

BigFloat::BigFloat(const BigInt& n) : mantissa_(n), exp_(0) {
    normalize();
}

BigFloat::BigFloat(const BigRat& r) {
    int p = g_precision;
    BigInt n = r.numerator().abs();
    BigInt d = r.denominator();
    BigInt scale = BigInt::pow(BigInt(10), BigInt(p + 10));
    BigInt m = n * scale / d;
    mantissa_ = m;
    exp_ = -(p + 10);
    if (r.numerator() < BigInt(0)) mantissa_ = -mantissa_;
    normalize();
}

BigFloat::BigFloat(const std::string& s) {
    if (s.empty() || s == "0") {
        mantissa_ = BigInt(0);
        exp_ = 0;
        return;
    }
    std::string clean;
    bool neg = false;
    int dot_pos = -1;
    int start = 0;
    if (s[0] == '-') { neg = true; start = 1; }
    else if (s[0] == '+') { start = 1; }
    int64_t sci_exp = 0;
    for (int i = start; i < (int)s.size(); i++) {
        if (s[i] == '.') { dot_pos = i - start; continue; }
        if (s[i] == 'e' || s[i] == 'E') {
            std::string exp_str = s.substr(i + 1);
            try { sci_exp = std::stoll(exp_str); } catch (...) { sci_exp = 0; }
            break;
        }
        if (s[i] >= '0' && s[i] <= '9') clean += s[i];
    }
    if (clean.empty()) {
        mantissa_ = BigInt(0);
        exp_ = 0;
        return;
    }
    mantissa_ = BigInt(clean);
    if (neg) mantissa_ = -mantissa_;
    int frac_digits = (dot_pos >= 0) ? (int)clean.size() - dot_pos : 0;
    exp_ = -frac_digits + sci_exp;
    normalize();
}

void BigFloat::normalize(int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (mantissa_.is_zero()) {
        exp_ = 0;
        return;
    }
    std::string s = mantissa_.abs().to_string();
    int digits = (int)s.size();
    if (digits > p + 5) {
        BigInt round_div = BigInt::pow(BigInt(10), BigInt(digits - p - 2));
        mantissa_ = mantissa_ / round_div;
        exp_ += digits - p - 2;
    }
    s = mantissa_.abs().to_string();
    digits = (int)s.size();
    while (digits > 1 && s.back() == '0') {
        mantissa_ = mantissa_ / BigInt(10);
        exp_ += 1;
        s = mantissa_.abs().to_string();
        digits = (int)s.size();
    }
}

BigFloat BigFloat::mul_int(const BigInt& n) const {
    BigFloat r;
    r.mantissa_ = mantissa_ * n;
    r.exp_ = exp_;
    r.normalize();
    return r;
}

BigFloat BigFloat::operator+(const BigFloat& o) const {
    if (is_zero()) return o;
    if (o.is_zero()) return *this;
    BigFloat a = *this, b = o;
    int64_t diff = a.exp_ - b.exp_;
    if (diff > 0) {
        a.mantissa_ = a.mantissa_ * BigInt::pow(BigInt(10), BigInt(diff));
        a.exp_ = b.exp_;
    } else if (diff < 0) {
        b.mantissa_ = b.mantissa_ * BigInt::pow(BigInt(10), BigInt(-diff));
        b.exp_ = a.exp_;
    }
    BigFloat r;
    r.mantissa_ = a.mantissa_ + b.mantissa_;
    r.exp_ = a.exp_;
    r.normalize();
    return r;
}

BigFloat BigFloat::operator-(const BigFloat& o) const {
    return *this + (-o);
}

BigFloat BigFloat::operator*(const BigFloat& o) const {
    BigFloat r;
    r.mantissa_ = mantissa_ * o.mantissa_;
    r.exp_ = exp_ + o.exp_;
    r.normalize();
    return r;
}

BigFloat BigFloat::operator/(const BigFloat& o) const {
    if (o.is_zero()) throw CalcError("Division by zero");
    if (is_zero()) return BigFloat(0);
    int p = g_precision;
    BigInt a = mantissa_.abs();
    BigInt b = o.mantissa_.abs();
    std::string sa = a.to_string();
    std::string sb = b.to_string();
    int extra = p + 10 - (int)sa.size() + (int)sb.size();
    if (extra < 0) extra = 0;
    BigInt scale = BigInt::pow(BigInt(10), BigInt(extra));
    BigInt q = a * scale / b;
    BigFloat r;
    r.mantissa_ = q;
    if (mantissa_.negative() != o.mantissa_.negative()) r.mantissa_ = -r.mantissa_;
    r.exp_ = exp_ - o.exp_ - extra;
    r.normalize();
    return r;
}

BigFloat BigFloat::operator-() const {
    BigFloat r = *this;
    r.mantissa_ = -r.mantissa_;
    return r;
}

bool BigFloat::operator==(const BigFloat& o) const {
    if (is_zero() && o.is_zero()) return true;
    BigFloat diff = *this - o;
    return diff.is_zero();
}

bool BigFloat::operator<(const BigFloat& o) const {
    BigFloat diff = *this - o;
    return diff.sign() < 0;
}

BigFloat BigFloat::abs() const {
    BigFloat r = *this;
    r.mantissa_ = r.mantissa_.abs();
    return r;
}

BigFloat BigFloat::round() const {
    BigFloat half("0.5");
    if (mantissa_.is_zero()) return BigFloat(0);
    if (*this >= BigFloat(0)) {
        std::string s = (*this + half).to_string(0);
        size_t dot = s.find('.');
        if (dot != std::string::npos) s = s.substr(0, dot);
        return BigFloat(s);
    } else {
        std::string s = (*this - half).to_string(0);
        size_t dot = s.find('.');
        if (dot != std::string::npos) s = s.substr(0, dot);
        return BigFloat(s);
    }
}

int BigFloat::sign() const {
    return mantissa_.sign();
}

bool BigFloat::is_zero() const {
    return mantissa_.is_zero();
}

std::string BigFloat::to_string(int prec) const {
    int p = prec > 0 ? prec : g_precision;
    if (is_zero()) return "0";
    std::string s = mantissa_.abs().to_string();
    bool neg = mantissa_.negative();
    int64_t point_pos = (int64_t)s.size() + exp_;
    if (point_pos <= 0) {
        s = std::string(-point_pos, '0') + s;
        point_pos = 0;
    }
    if ((int64_t)s.size() > p + 5) {
        s = s.substr(0, p + 5);
    }
    while ((int)s.size() > 1 && s.back() == '0') s.pop_back();
    std::string result;
    if (neg) result += "-";
    if (point_pos >= (int64_t)s.size()) {
        result += s + std::string(point_pos - s.size(), '0');
    } else if (point_pos <= 0) {
        result += "0." + std::string(-point_pos, '0') + s;
    } else {
        result += s.substr(0, point_pos);
        if (point_pos < (int64_t)s.size()) {
            std::string frac = s.substr(point_pos);
            while (!frac.empty() && frac.back() == '0') frac.pop_back();
            if (!frac.empty()) result += "." + frac;
        }
    }
    return result;
}

BigFloat BigFloat::pi(int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (!pi_cache_.is_zero()) {
        std::string cached = pi_cache_.mantissa_.abs().to_string();
        if ((int)cached.size() >= p + 5) return pi_cache_;
    }
    BigFloat result("3.14159265358979323846264338327950288419716939937510"
                    "58209749445923078164062862089986280348253421170679"
                    "82148086513282306647093844609550582231725359408128"
                    "48111745028410270193852110555964462294895493038196"
                    "44288109756659334461284756482337867831652712019091"
                    "45648566923460348610454326648213393607260249141273"
                    "72458700660631558817488152092096282925409171536436"
                    "78925903600113305305488204665213841469519415116094"
                    "33057270365759591953092186117381932611793105118548"
                    "07446237996274956735188575272489122793818301194912"
                    "98336733624406566430860213949463952247371907021798"
                    "60943702770539217176293176752384674818467669405132"
                    "00056812714526356082778577134275778960917363717872"
                    "14684409012249534301465495853710507922796892589235"
                    "42019956112129021960864034418159813629774771309960"
                    "51870721134999999837297804995105973173281609631859"
                    "50244594553469083026425223082533446850352619311881"
                    "71010003137838752886587533208381420617177669147303"
                    "59825349042875546873115956286388235378759375195778"
                    "18577805321712268066130019278766111959092164201989");
    result.normalize();
    pi_cache_ = result;
    return result;
}

BigFloat BigFloat::pi_cache_;

BigFloat BigFloat::e_val(int prec) {
    (void)prec;
    BigFloat result("2.71828182845904523536028747135266249775724709369995"
                    "95749669676277240766303535475945713821785251664274"
                    "27466391932003059921817413596629043572900334295260"
                    "59563073813232862794349076323382988075319525101901"
                    "15738341879307021540891499348841675092447614606680"
                    "82264800168447754436879018635801632420137296269932"
                    "21306085633700065401730071141887629930043714632519"
                    "61044525261168854124305860252864764155200402149057"
                    "51691526854844844671501682405193704851668061711716"
                    "40608671702026168831385833057604253664891762158018"
                    "87432723096684921920480252191079553024947409943897"
                    "68446044663646995812867177816731662014339705693669"
                    "60751394891622636269704454045870861781086859594660"
                    "49989576838164179056504656595076863944246566528968"
                    "83400690325860279333967549850449561127943440660413"
                    "27647665254953953747802849184931854144540363882241"
                    "09030774254654883490674188492963020457794461058658"
                    "41986830789046298217819489087242103068951074264884"
                    "42483057660262486267917342071009440583248632658816"
                    "82444179450292591792500474184264826044538878670842");
    result.normalize();
    return result;
}

BigFloat BigFloat::ln2(int prec) {
    (void)prec;
    BigFloat result("0.69314718055994530941723212145817656807550013436025"
                    "52541206800094933936219696947156058633269964186875420"
                    "01481021057084131167231602349227453264004854979355350"
                    "72647004973871676839634264483992534866057772662900081"
                    "68641668842476875695814061068073374857106352946505688"
                    "41835878486867987802455653773465564473277155004028677"
                    "56871093489642900784691064063669682871926170824016758"
                    "07205292628683441031386863263886657035978547197529474"
                    "54692654917351096305584763518543696792082993785031866"
                    "39063025950296520483030942949726192887101164284733197"
                    "61031463977426364404276454030784052931462001232986877"
                    "62841552891023094674487696857562166564025143991146677"
                    "02153023579257864834283958839754066539895893591643710"
                    "63967663987191009981493200696879054094984292202878568"
                    "31284229103462570768749673289687074882043862766564327"
                    "94001640602559602774353774050654467704532479815983865"
                    "61062394967241550376453077196654064230491064539286878"
                    "40401965661693853026437902825239079194593944069061270"
                    "84839957897102089546546529694602084029963853958063717"
                    "64479940301445600035663404064991597687364408753349967");
    result.normalize();
    return result;
}

BigFloat BigFloat::exp(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (x.is_zero()) return BigFloat(1);
    BigFloat ln2_val = BigFloat::ln2(p + 20);
    int64_t k = 0;
    if (!ln2_val.is_zero()) {
        BigFloat kf = x / ln2_val;
        std::string ks = kf.to_string(20);
        size_t dot = ks.find('.');
        if (dot != std::string::npos) ks = ks.substr(0, dot);
        if (!ks.empty()) {
            try { k = std::stoll(ks); } catch (...) { k = 0; }
        }
    }
    BigFloat r = x - ln2_val * BigFloat(k);
    while (r > ln2_val / BigFloat(2)) { r = r - ln2_val; k++; }
    while (r < -ln2_val / BigFloat(2)) { r = r + ln2_val; k--; }
    BigFloat result(1);
    BigFloat term(1);
    for (int64_t i = 1; i <= p + 100; i++) {
        term = term * r / BigFloat(i);
        result = result + term;
        if (term.abs() < BigFloat("1e-" + std::to_string(p + 10))) break;
    }
    for (int64_t i = 0; i < k; i++) result = result * BigFloat(2);
    for (int64_t i = 0; i < -k; i++) result = result / BigFloat(2);
    result.normalize();
    return result;
}

BigFloat BigFloat::ln(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (x.is_zero()) throw CalcError("ln(0) is undefined");
    if (x < BigFloat(0)) throw CalcError("ln of negative number");
    if (x == BigFloat(1)) return BigFloat(0);
    int64_t k = 0;
    BigFloat m = x;
    BigFloat half_sqrt2 = BigFloat("0.70710678118654752440084436210484903928483593768847"
                                    "40365883398689953662392310535194251937671638207863");
    while (m > BigFloat(2)) { m = m / BigFloat(2); k++; }
    while (m < half_sqrt2) { m = m * BigFloat(2); k--; }
    BigFloat y = (m - BigFloat(1)) / (m + BigFloat(1));
    BigFloat y2 = y * y;
    BigFloat term = y;
    BigFloat sum = y;
    for (int64_t i = 1; i <= p + 100; i++) {
        term = term * y2;
        BigFloat contrib = term / BigFloat(2 * i + 1);
        sum = sum + contrib;
        if (contrib.abs() < BigFloat("1e-" + std::to_string(p + 10))) break;
    }
    BigFloat ln2_val = BigFloat::ln2(p + 20);
    return sum * BigFloat(2) + ln2_val * BigFloat(k);
}

BigFloat BigFloat::log10(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    return ln(x, p + 5) / ln(BigFloat(10), p + 5);
}

BigFloat BigFloat::log2(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    return ln(x, p + 5) / ln2(p + 5);
}

BigFloat BigFloat::log_base(const BigFloat& x, const BigFloat& base, int prec) {
    int p = prec > 0 ? prec : g_precision;
    return ln(x, p + 5) / ln(base, p + 5);
}

BigFloat BigFloat::sin_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (x.is_zero()) return BigFloat(0);
    BigFloat pi_val = BigFloat::pi(p + 20);
    BigFloat halfpi = pi_val / BigFloat(2);
    BigFloat twopi = pi_val * BigFloat(2);
    BigFloat r = x;
    while (r > twopi) r = r - twopi;
    while (r < BigFloat(0)) r = r + twopi;
    int quadrant = 0;
    while (r > halfpi) {
        r = r - halfpi;
        quadrant++;
    }
    if (quadrant % 2 == 1) {
        r = halfpi - r;
    }
    BigFloat result = r;
    BigFloat term = r;
    BigFloat r2 = r * r;
    for (int64_t i = 1; i <= p + 100; i++) {
        term = term * r2 * BigFloat(-1) / BigFloat((2 * i) * (2 * i + 1));
        result = result + term;
        if (term.abs() < BigFloat("1e-" + std::to_string(p + 10))) break;
    }
    if (quadrant % 4 >= 2) result = -result;
    result.normalize();
    return result;
}

BigFloat BigFloat::cos_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (x.is_zero()) return BigFloat(1);
    BigFloat pi_val = BigFloat::pi(p + 20);
    BigFloat halfpi = pi_val / BigFloat(2);
    BigFloat twopi = pi_val * BigFloat(2);
    BigFloat r = x;
    while (r > twopi) r = r - twopi;
    while (r < BigFloat(0)) r = r + twopi;
    int quadrant = 0;
    while (r > halfpi) {
        r = r - halfpi;
        quadrant++;
    }
    if (quadrant % 2 == 1) {
        r = halfpi - r;
    }
    BigFloat result(1);
    BigFloat term(1);
    BigFloat r2 = r * r;
    for (int64_t i = 1; i <= p + 100; i++) {
        term = term * r2 * BigFloat(-1) / BigFloat((2 * i - 1) * (2 * i));
        result = result + term;
        if (term.abs() < BigFloat("1e-" + std::to_string(p + 10))) break;
    }
    if (quadrant % 4 >= 2) result = -result;
    result.normalize();
    return result;
}

BigFloat BigFloat::tan_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    BigFloat s = sin_val(x, p + 5);
    BigFloat c = cos_val(x, p + 5);
    if (c.is_zero()) throw CalcError("tan is undefined at this point");
    return s / c;
}

BigFloat BigFloat::asin_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (x > BigFloat(1) || x < BigFloat(-1))
        throw CalcError("asin argument out of range [-1, 1]");
    if (x.is_zero()) return BigFloat(0);
    if (x == BigFloat(1)) return BigFloat::pi(p) / BigFloat(2);
    if (x == BigFloat(-1)) return -BigFloat::pi(p) / BigFloat(2);
    BigFloat abs_x = x.abs();
    if (abs_x < BigFloat("0.5")) {
        BigFloat result = x;
        BigFloat term = x;
        BigFloat x2 = x * x;
        for (int64_t i = 1; i <= p + 100; i++) {
            term = term * x2 * BigFloat(2 * i - 1) * BigFloat(2 * i - 1) /
                   (BigFloat(2 * i) * BigFloat(2 * i + 1));
            result = result + term;
            if (term.abs() < BigFloat("1e-" + std::to_string(p + 10))) break;
        }
        return result;
    }
    BigFloat pi_val = BigFloat::pi(p + 20);
    BigFloat half = BigFloat("0.5");
    BigFloat one = BigFloat(1);
    BigFloat y = BigFloat::sqrt_val(one - x * x, p + 10);
    int sign = x.sign();
    return pi_val / BigFloat(2) - atan_val(y, p + 10) * BigFloat(sign);
}

BigFloat BigFloat::acos_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (x > BigFloat(1) || x < BigFloat(-1))
        throw CalcError("acos argument out of range [-1, 1]");
    BigFloat pi_val = BigFloat::pi(p + 10);
    return pi_val / BigFloat(2) - asin_val(x, p + 10);
}

BigFloat BigFloat::atan_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (x.is_zero()) return BigFloat(0);
    BigFloat pi_val = BigFloat::pi(p + 20);
    bool negate = false;
    BigFloat a = x;
    if (a < BigFloat(0)) { a = -a; negate = true; }
    bool invert = false;
    if (a > BigFloat(1)) { a = BigFloat(1) / a; invert = true; }
    BigFloat result = a;
    BigFloat term = a;
    BigFloat a2 = a * a;
    for (int64_t i = 1; i <= p + 200; i++) {
        term = term * a2 * BigFloat(-1);
        BigFloat contrib = term / BigFloat(2 * i + 1);
        result = result + contrib;
        if (contrib.abs() < BigFloat("1e-" + std::to_string(p + 10))) break;
    }
    if (invert) result = pi_val / BigFloat(2) - result;
    if (negate) result = -result;
    result.normalize();
    return result;
}

BigFloat BigFloat::sinh_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    BigFloat ep = exp(x, p + 10);
    BigFloat em = exp(-x, p + 10);
    return (ep - em) / BigFloat(2);
}

BigFloat BigFloat::cosh_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    BigFloat ep = exp(x, p + 10);
    BigFloat em = exp(-x, p + 10);
    return (ep + em) / BigFloat(2);
}

BigFloat BigFloat::tanh_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    BigFloat s = sinh_val(x, p + 5);
    BigFloat c = cosh_val(x, p + 5);
    return s / c;
}

BigFloat BigFloat::sqrt_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (x.is_zero()) return BigFloat(0);
    if (x < BigFloat(0)) throw CalcError("sqrt of negative number");
    double approx = 1.0;
    std::string s = x.to_string(20);
    try { approx = std::stod(s); } catch (...) { approx = 1.0; }
    if (approx <= 0) approx = 1.0;
    approx = std::sqrt(approx);
    BigFloat result(approx);
    result.normalize();
    for (int i = 0; i < p + 50; i++) {
        BigFloat next = (result + x / result) / BigFloat(2);
        BigFloat diff = (next - result).abs();
        result = next;
        if (diff.is_zero()) break;
        std::string ds = diff.to_string();
        if ((int)ds.size() < 2 && ds[0] == '0') break;
    }
    result.normalize();
    return result;
}

BigFloat BigFloat::cbrt_val(const BigFloat& x, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (x.is_zero()) return BigFloat(0);
    bool neg = x < BigFloat(0);
    BigFloat a = neg ? -x : x;
    double approx = 1.0;
    std::string s = a.to_string(20);
    try { approx = std::stod(s); } catch (...) { approx = 1.0; }
    if (approx <= 0) approx = 1.0;
    approx = std::cbrt(approx);
    BigFloat result(approx);
    result.normalize();
    for (int i = 0; i < p + 50; i++) {
        BigFloat r2 = result * result;
        BigFloat next = (BigFloat(2) * result + a / r2) / BigFloat(3);
        BigFloat diff = (next - result).abs();
        result = next;
        if (diff.is_zero()) break;
    }
    if (neg) result = -result;
    result.normalize();
    return result;
}

BigFloat BigFloat::nrt_val(const BigFloat& x, int64_t n, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (n <= 0) throw CalcError("n-th root requires positive n");
    if (n == 1) return x;
    if (n == 2) return sqrt_val(x, p);
    if (n == 3) return cbrt_val(x, p);
    if (x.is_zero()) return BigFloat(0);
    bool neg = (n % 2 == 1) && x < BigFloat(0);
    BigFloat a = neg ? -x : x;
    if (x < BigFloat(0) && n % 2 == 0)
        throw CalcError("even root of negative number");
    double approx = 1.0;
    std::string s = a.to_string(20);
    try { approx = std::stod(s); } catch (...) { approx = 1.0; }
    if (approx <= 0) approx = 1.0;
    approx = std::pow(approx, 1.0 / n);
    BigFloat result(approx);
    result.normalize();
    BigFloat nf(n);
    BigFloat nm1(n - 1);
    for (int i = 0; i < p + 50; i++) {
        BigFloat rn = pow_val(result, nf, p + 10);
        BigFloat next = (nm1 * result + a / rn) / nf;
        BigFloat diff = (next - result).abs();
        result = next;
        if (diff.is_zero()) break;
    }
    if (neg) result = -result;
    result.normalize();
    return result;
}

BigFloat BigFloat::pow_val(const BigFloat& base, const BigFloat& exponent, int prec) {
    int p = prec > 0 ? prec : g_precision;
    if (exponent.is_zero()) return BigFloat(1);
    if (base.is_zero()) {
        if (exponent > BigFloat(0)) return BigFloat(0);
        throw CalcError("0 raised to negative power");
    }
    if (base < BigFloat(0)) {
        BigFloat exp_rational = exponent;
        std::string es = exponent.to_string(30);
        size_t dot = es.find('.');
        if (dot != std::string::npos) {
            bool all_zero = true;
            for (size_t i = dot + 1; i < es.size(); i++) {
                if (es[i] != '0') { all_zero = false; break; }
            }
            if (all_zero) {
                int64_t n = 0;
                try { n = std::stoll(es.substr(0, dot)); } catch (...) {}
                if (n % 2 == 0) return exp(exponent * ln(-base, p + 10), p + 10);
                return -exp(exponent * ln(-base, p + 10), p + 10);
            }
        }
        throw CalcError("Cannot compute non-integer power of negative number (use complex)");
    }
    BigFloat ln_base = ln(base, p + 10);
    return exp(exponent * ln_base, p + 10);
}

BigFloat BigFloat::factorial(int64_t n) {
    if (n < 0) throw CalcError("Factorial of negative number");
    BigInt result = BigInt::factorial(n);
    return BigFloat(result);
}
