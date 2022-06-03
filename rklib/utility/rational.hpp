#ifndef RK_RATIONAL_HPP
#define RK_RATIONAL_HPP

#include <cassert>
#include <iostream>

namespace rklib {

namespace rational_internal {

template <class T>
T gcd(T a, T b) {
    if (b == T(0)) {
        if (a < T(0)) a = -a;
        return a;
    }
    return gcd(b, a % b);
}

}  // namespace rational_internal

template <class T = long long>
struct Rational {
   public:
    Rational() : num(0), den(1){};
    Rational(T n) : num(n), den(1) {}
    Rational(T n, T d) : num(n), den(d) {
        assert(den != T(0));
        if (den < 0) num = -num, den = -den;
        T g = rational_internal::gcd(num, den);
        num /= g;
        den /= g;
    }

    T numerator() { return num; }
    T denominator() { return den; }

    Rational& operator+=(const Rational& rhs) {
        T r_den = rhs.den;
        T g = rational_internal::gcd(den, r_den);
        den /= g;
        num = num * (r_den / g) + rhs.num * den;
        g = rational_internal::gcd(num, g);
        num /= g;
        den *= r_den / g;
        return *this;
    }
    Rational& operator-=(const Rational& rhs) {
        T r_den = rhs.den;
        T g = rational_internal::gcd(den, r_den);
        den /= g;
        num = num * (r_den / g) - rhs.num * den;
        g = rational_internal::gcd(num, g);
        num /= g;
        den *= r_den / g;
        return *this;
    }
    Rational& operator*=(const Rational& rhs) {
        T r_num = rhs.num;
        T r_den = rhs.den;

        T g1 = rational_internal::gcd(num, r_den);
        T g2 = rational_internal::gcd(r_num, den);
        num = (num / g1) * (r_num / g2);
        den = (den / g2) * (r_den / g1);
        return *this;
    }
    Rational& operator/=(const Rational& rhs) {
        T r_num = rhs.num;
        T r_den = rhs.den;

        T zero(0);

        assert(r_num != zero);
        if (num == zero) return *this;

        T gcd1 = rational_internal::gcd(num, r_num);
        T gcd2 = rational_internal::gcd(r_den, den);
        num = (num / gcd1) * (r_den / gcd2);
        den = (den / gcd2) * (r_num / gcd1);

        if (den < zero) {
            num = -num;
            den = -den;
        }
        return *this;
    }

    friend Rational operator+(const Rational& lhs, const Rational& rhs) {
        return Rational(lhs) += rhs;
    }
    friend Rational operator-(const Rational& lhs, const Rational& rhs) {
        return Rational(lhs) -= rhs;
    }
    friend Rational operator*(const Rational& lhs, const Rational& rhs) {
        return Rational(lhs) *= rhs;
    }
    friend Rational operator/(const Rational& lhs, const Rational& rhs) {
        return Rational(lhs) /= rhs;
    }

    Rational operator+() const { return *this; }
    Rational operator-() const { return Rational{-num, den}; }

    Rational& operator++() {
        num += den;
        return *this;
    }
    Rational& operator--() {
        num -= den;
        return *this;
    }

    friend bool operator==(const Rational& lhs, const Rational& rhs) {
        return lhs.num == rhs.num && lhs.den == rhs.den;
    }
    friend bool operator!=(const Rational& lhs, const Rational& rhs) {
        return !(lhs == rhs);
    }
    friend bool operator<(const Rational& lhs, const Rational& rhs) {
        return lhs.num * rhs.den < rhs.num * lhs.den;
    }
    friend bool operator>(const Rational& lhs, const Rational& rhs) {
        return rhs < lhs;
    }
    friend bool operator<=(const Rational& lhs, const Rational& rhs) {
        return !(lhs > rhs);
    }
    friend bool operator>=(const Rational& lhs, const Rational& rhs) {
        return rhs <= lhs;
    }

    friend std::ostream& operator<<(std::ostream& os, const Rational& r) {
        if (r.den == T(1)) return os << r.num;
        return os << r.num << "/" << r.den;
    }

    operator double() const { return double(num) / double(den); }
    operator long double() const { return (long double)num / (long double)den; }

   private:
    T num, den;
};

}  // namespace rklib

#endif  // RK_RATIONAL_HPP