#ifndef RK_BIG_INTEGER_HPP
#define RK_BIG_INTEGER_HPP

#include <iostream>
#include <rklib/utility/utility.hpp>
#include <string>
#include <type_traits>
#include <vector>

namespace rklib {

struct BigInt {
   public:
    BigInt() : BigInt(0) {}
    template <typename T, std::enable_if_t<std::is_integral<T>::value,
                                           std::nullptr_t> = nullptr>
    BigInt(T x) : BigInt(std::vector<T>(1, x)) {}
    template <typename T, std::enable_if_t<std::is_integral<T>::value,
                                           std::nullptr_t> = nullptr>
    BigInt(std::vector<T> a) {
        while (a.size() >= 2 && a.back() == 0) a.pop_back();
        if (a.size() == 1 && a[0] == 0) {
            sgn = 0;
            v = {0};
            return;
        }
        sgn = (a.back() > 0 ? 1 : -1);
        if (sgn < 0) {
            std::transform(a.begin(), a.end(), a.begin(),
                           [](T x) { return -x; });
        }
        v = normalize(a);
    }
    BigInt(std::string s) : sgn(1) {
        if (s[0] == '-') {
            s.erase(s.begin());
            sgn = -1;
        }
        if (s == "0") {
            sgn = 0;
            v = {0};
            return;
        }
        for (int r = int(s.size()), l = std::max(0, r - log10_base); r > 0;
             r = l, l = std::max(0, r - log10_base)) {
            int tmp = 0;
            for (int i = l; i < r; ++i) {
                tmp = tmp * 10 + (s[i] - '0');
            }
            v.push_back(tmp);
        }
    }

    BigInt& operator+=(const BigInt& rhs) {
        int r_sgn = rhs.sgn;
        if (r_sgn == 0) return *this;
        if (sgn == 0) return *this = rhs;
        if (v.size() < rhs.v.size()) v.resize(rhs.v.size(), 0);
        for (size_t i = 0; i < rhs.v.size(); i++) {
            if (sgn == r_sgn)
                v[i] += rhs.v[i];
            else
                v[i] -= rhs.v[i];
        }
        normalize();
        return *this;
    }
    BigInt& operator-=(const BigInt& rhs) {
        int r_sgn = rhs.sgn;
        if (r_sgn == 0) return *this;
        if (sgn == 0) return *this = -rhs;
        if (v.size() < rhs.v.size()) v.resize(rhs.v.size(), 0);
        for (size_t i = 0; i < rhs.v.size(); i++) {
            if (sgn == r_sgn)
                v[i] -= rhs.v[i];
            else
                v[i] += rhs.v[i];
        }
        normalize();
        return *this;
    }
    BigInt& operator*=(const BigInt& rhs) {
        int r_sgn = rhs.sgn;
        if (sgn == 0) return *this;
        if (r_sgn == 0) return *this = rhs;

        int res_sgn = sgn * r_sgn;

        if (int(std::min(v.size(), rhs.v.size())) <= naive_mul_threshold) {
            v = multiply_naive(v, rhs.v);
        } else {
            std::vector<long long> a(v.size()), b(rhs.v.size());
            for (size_t i = 0; i < a.size(); i++) {
                a[i] = (long long)v[i];
            }
            for (size_t i = 0; i < b.size(); i++) {
                b[i] = (long long)rhs.v[i];
            }
            int sz = 1;
            while (sz < int(a.size()) || sz < int(b.size())) sz *= 2;
            a.resize(sz, 0);
            b.resize(sz, 0);

            this->v = normalize(multiply_karatsuba(a, b));
        }

        this->sgn = res_sgn;

        return *this;
    }
    BigInt& operator/=(const BigInt& rhs) {
        int r_sgn = rhs.sgn;
        assert(r_sgn != 0);
        if (sgn == 0) return *this;
        if (rhs == BigInt(1) || rhs == BigInt(-1)) {
            sgn *= r_sgn;
            return *this;
        }

        int res_sgn = sgn * r_sgn;

        (*this) = divide_newton_raphson(*this, rhs);

        this->sgn = (this->v.back() == 0 ? 0 : res_sgn);

        return *this;
    }
    BigInt& operator%=(const BigInt& rhs) {
        assert(rhs.sgn != 0);

        return *this = (*this) - (*this) / rhs * rhs;
    }

    friend BigInt operator+(const BigInt& lhs, const BigInt& rhs) {
        return BigInt(lhs) += rhs;
    }
    friend BigInt operator-(const BigInt& lhs, const BigInt& rhs) {
        return BigInt(lhs) -= rhs;
    }
    friend BigInt operator*(const BigInt& lhs, const BigInt& rhs) {
        return BigInt(lhs) *= rhs;
    }
    friend BigInt operator/(const BigInt& lhs, const BigInt& rhs) {
        return BigInt(lhs) /= rhs;
    }
    friend BigInt operator%(const BigInt& lhs, const BigInt& rhs) {
        return BigInt(lhs) %= rhs;
    }

    BigInt operator+() const { return *this; }
    BigInt operator-() const { return {-sgn, v}; }

    BigInt& operator++() { return *this += 1; }
    BigInt& operator--() { return *this -= 1; }

    friend bool operator==(const BigInt& lhs, const BigInt& rhs) {
        if (lhs.sgn != rhs.sgn || lhs.v.size() != rhs.v.size()) return false;
        for (size_t i = 0; i < lhs.v.size(); i++) {
            if (lhs.v[i] != rhs.v[i]) return false;
        }
        return true;
    }
    friend bool operator!=(const BigInt& lhs, const BigInt& rhs) {
        return !(lhs == rhs);
    }
    friend bool operator<(const BigInt& lhs, const BigInt& rhs) {
        int l_sgn = lhs.sgn, r_sgn = rhs.sgn;
        if (l_sgn < r_sgn) return true;
        if (l_sgn > r_sgn) return false;

        int nl = lhs.v.size(), nr = rhs.v.size();
        if (l_sgn * nl < r_sgn * nr) return true;
        if (l_sgn * nl > r_sgn * nr) return false;

        for (int i = nl - 1; i >= 0; i--) {
            if (l_sgn * lhs.v[i] < r_sgn * rhs.v[i]) return true;
            if (l_sgn * lhs.v[i] > r_sgn * rhs.v[i]) return false;
        }

        return false;
    }
    friend bool operator>(const BigInt& lhs, const BigInt& rhs) {
        return rhs < lhs;
    }
    friend bool operator<=(const BigInt& lhs, const BigInt& rhs) {
        return !(lhs > rhs);
    }
    friend bool operator>=(const BigInt& lhs, const BigInt& rhs) {
        return rhs <= lhs;
    }

    friend std::istream& operator>>(std::istream& is, BigInt& rhs) {
        std::string s;
        is >> s;
        rhs = BigInt(s);
        return is;
    }
    friend std::ostream& operator<<(std::ostream& os, const BigInt& rhs) {
        if (rhs.sgn < 0) os << "-";
        for (int i = int(rhs.v.size()) - 1; i >= 0; i--) {
            if (i == int(rhs.v.size()) - 1) {
                os << rhs.v[i];
            } else {
                os << std::to_string(rhs.v[i] + base).substr(1, log10_base);
            }
        }
        return os;
    }

    operator double() const {
        double res = 0;
        for (int i = v.size() - 1; i >= 0; i--) {
            res = res * base + v[i];
        }
        res *= sgn;
        return res;
    }
    operator long double() const {
        long double res = 0;
        for (int i = v.size() - 1; i >= 0; i--) {
            res = res * base + v[i];
        }
        res *= sgn;
        return res;
    }

   private:
    static constexpr int base = 10000, log10_base = 4;
    static constexpr int naive_mul_threshold = 75;
    int sgn;
    std::vector<int> v;

    BigInt(int sgn, std::vector<int> v) : sgn(sgn), v(v) {}

    void normalize() {
        while (v.size() >= 2 && v.back() == 0) v.pop_back();

        if (v.back() < 0) {
            sgn = -sgn;
            std::transform(v.begin(), v.end(), v.begin(),
                           [](int x) { return -x; });
        }

        int carry = 0;
        for (size_t i = 0; i < v.size(); i++) {
            v[i] += carry;
            if (0 <= v[i] && v[i] < base) {
                carry = 0;
                continue;
            }
            carry = div_floor(v[i], base);
            v[i] -= carry * base;
        }
        while (carry > 0) {
            v.push_back(carry % base);
            carry /= base;
        }

        while (v.size() >= 2 && v.back() == 0) v.pop_back();
        if (v.size() == 1 && v[0] == 0) sgn = 0;
    }
    template <class T>
    static std::vector<int> normalize(const std::vector<T>& c) {
        std::vector<int> res(c.size());
        T carry = 0;
        for (size_t i = 0; i < c.size(); i++) {
            T tmp = c[i];
            tmp += carry;
            if (0 <= tmp && tmp < T(base)) {
                carry = 0;
                res[i] = int(tmp);
                continue;
            }
            carry = div_floor(tmp, T(base));
            res[i] = int(tmp - carry * base);
        }
        while (carry > 0) {
            res.push_back(int(carry % base));
            carry /= base;
        }

        while (res.size() >= 2 && res.back() == 0) res.pop_back();

        return res;
    }

    static std::vector<int> multiply_naive(const std::vector<int>& a,
                                           const std::vector<int>& b) {
        std::vector<long long> c(a.size() + b.size() - 1, 0);
        for (size_t i = 0; i < a.size(); i++) {
            for (size_t j = 0; j < b.size(); j++) {
                c[i + j] += a[i] * b[j];
            }
        }

        return normalize(c);
    }
    static std::vector<long long> multiply_karatsuba(std::vector<long long> a,
                                                     std::vector<long long> b) {
        int na = a.size(), nb = b.size();
        if (std::min(na, nb) <= naive_mul_threshold) {
            std::vector<long long> res(a.size() + b.size() - 1, 0);
            for (size_t i = 0; i < a.size(); i++) {
                for (size_t j = 0; j < b.size(); j++) {
                    res[i + j] += a[i] * b[j];
                }
            }
            return res;
        }

        int n = std::max(na, nb);
        if (na < n) a.resize(n, 0);
        if (nb < n) b.resize(n, 0);

        int k = n / 2;
        std::vector<long long> x(
            std::vector<long long>(a.begin() + k, a.end()));
        std::vector<long long> y(
            std::vector<long long>(a.begin(), a.begin() + k));
        std::vector<long long> z(
            std::vector<long long>(b.begin() + k, b.end()));
        std::vector<long long> w(
            std::vector<long long>(b.begin(), b.begin() + k));

        std::vector<long long> xz = multiply_karatsuba(x, z);
        std::vector<long long> yw = multiply_karatsuba(y, w);

        for (int i = 0; i < k; i++) {
            x[i] += y[i];
            z[i] += w[i];
        }
        std::vector<long long> t = multiply_karatsuba(x, z);
        for (size_t i = 0; i < xz.size(); i++) {
            t[i] -= xz[i] + yw[i];
        }

        a = std::vector<long long>(2 * k, 0);
        a.insert(a.end(), xz.begin(), xz.end());
        b = std::vector<long long>(k, 0);
        b.insert(b.end(), t.begin(), t.end());

        for (size_t i = 0; i < b.size(); i++) {
            a[i] += b[i];
        }
        for (size_t i = 0; i < yw.size(); i++) {
            a[i] += yw[i];
        }

        return a;
    }

    static void shift_r(BigInt& a, size_t len) {
        if (a == BigInt(0)) return;
        if (a.v.size() <= len) {
            a = 0;
            return;
        }
        a.v.erase(a.v.begin(), a.v.begin() + len);
    }

    static BigInt divide_newton_raphson(BigInt a, BigInt b) {
        a.sgn = b.sgn = 1;
        if (a < b) return 0;
        size_t len = a.v.size() + 2;
        BigInt x = base * base, prv = x;
        while (true) {
            BigInt tmp = x * x * b;
            shift_r(tmp, len);
            x += x;
            x -= tmp;
            if (x.v.size() > len) shift_r(x, x.v.size() - len);

            if ((x - prv).v.size() <= 1) break;
            prv = x;
        }
        BigInt res = a * x;
        res.v.erase(res.v.begin(), res.v.begin() + len);
        if (res * b + b <= a) ++res;

        return res;
    }
};

}  // namespace rklib

#endif  // RK_BIG_INTEGER_HPP