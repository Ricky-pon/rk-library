#ifndef RK_FORMAL_POWER_SERIES_HPP
#define RK_FORMAL_POWER_SERIES_HPP

#include <atcoder/convolution>
#include <rklib/math/factorial.hpp>
#include <rklib/utility/utility.hpp>
#include <vector>

namespace rklib {

template <class T>
struct FormalPowerSeries : public std::vector<T> {
   public:
    using fps = FormalPowerSeries;
    using std::vector<T>::vector;
    using std::vector<T>::operator=;

    fps& operator+=(const fps& rhs) {
        if (this->size() < rhs.size()) this->resize(rhs.size(), 0);
        for (size_t i = 0; i < this->size() && i < rhs.size(); i++) {
            (*this)[i] += rhs[i];
        }
        return *this;
    }
    fps& operator-=(const fps& rhs) {
        if (this->size() < rhs.size()) this->resize(rhs.size(), 0);
        for (size_t i = 0; i < this->size() && i < rhs.size(); i++) {
            (*this)[i] -= rhs[i];
        }
        return *this;
    }
    fps& operator*=(const fps& rhs) {
        *this = atcoder::convolution((*this), rhs);
        return *this;
    }
    fps& operator*=(const T& rhs) {
        std::transform((*this).begin(), (*this).end(), (*this).begin(),
                       [&](T x) { return x * rhs; });
        return *this;
    }
    fps& operator/=(const T& rhs) {
        T rhs_inv = T(1) / rhs;
        std::transform((*this).begin(), (*this).end(), (*this).begin(),
                       [&](T x) { return x * rhs_inv; });
        return *this;
    }
    fps& operator/=(const fps& rhs) {
        *this *= rhs.inv();
        this->resize(rhs.size());
        return *this;
    }

    fps operator+() const { return *this; }
    fps operator-() const {
        return FormalPowerSeries<T>(this->size(), 0) - (*this);
    }

    friend fps operator+(const fps& lhs, const fps& rhs) {
        return fps(lhs) += rhs;
    }
    friend fps operator-(const fps& lhs, const fps& rhs) {
        return fps(lhs) -= rhs;
    }
    friend fps operator*(const fps& lhs, const fps& rhs) {
        return fps(lhs) *= rhs;
    }
    friend fps operator*(const fps& lhs, const T& rhs) {
        return fps(lhs) *= rhs;
    }
    friend fps operator*(const T& lhs, const fps& rhs) {
        return fps(rhs) *= lhs;
    }
    friend fps operator/(const fps& lhs, const T& rhs) {
        return fps(lhs) /= rhs;
    }

    fps inv() const {
        size_t n = this->size();
        fps res = {1 / this->front()};
        for (size_t m = 1; m < n; m <<= 1) {
            fps f(this->begin(), this->begin() + std::min(n, 2 * m)), g(res);

            f.resize(2 * m, 0);
            g.resize(2 * m, 0);
            atcoder::internal::butterfly(f);
            atcoder::internal::butterfly(g);
            for (size_t i = 0; i < f.size(); i++) {
                f[i] *= g[i];
            }
            atcoder::internal::butterfly_inv(f);
            f /= 2 * m;

            f.erase(f.begin(), f.begin() + m);
            f.resize(2 * m, 0);
            atcoder::internal::butterfly(f);
            for (size_t i = 0; i < f.size(); i++) {
                f[i] *= -g[i];
            }
            atcoder::internal::butterfly_inv(f);
            f /= 2 * m;

            res.insert(res.end(), f.begin(), f.begin() + m);
        }
        res.resize(n);
        return res;
    }

    fps derivative() const {
        fps res(this->begin() + 1, this->end());
        res.push_back(0);
        for (size_t i = 0; i < res.size(); i++) {
            res[i] *= (i + 1);
        }
        return res;
    }

    fps integral() const {
        fps res(this->begin(), this->end() - 1);
        res.insert(res.begin(), 0);
        for (size_t i = 1; i < res.size(); i++) {
            res[i] *= fc.inv(i);
        }
        return res;
    }

    fps log() const {
        fps f = this->derivative();
        f /= (*this);
        return f.integral();
    }

    fps exp() const {
        size_t n = this->size();
        fps f = {1}, g = {1};
        for (size_t m = 1; m < n; m <<= 1) {
            fps f_ft(f);
            f_ft.resize(2 * m, 0);
            atcoder::internal::butterfly(f_ft);

            g.resize(2 * m, 0);
            atcoder::internal::butterfly(g);
            for (size_t i = 0; i < 2 * m; i++) {
                g[i] = 2 * g[i] - f_ft[i] * g[i] * g[i];
            }
            atcoder::internal::butterfly_inv(g);
            g /= 2 * m;
            g.resize(m);

            fps q(this->begin(), this->begin() + std::min(m, n));
            q.resize(2 * m, 0);
            q = q.derivative();

            fps w;
            {
                fps g_ft(g);
                g_ft.resize(2 * m, 0);
                atcoder::internal::butterfly(g_ft);
                w = g_ft;
                atcoder::internal::butterfly(w);
                for (size_t i = 0; i < 2 * m; i++) {
                    w[i] = f_ft[i] * g_ft[i];
                }
                atcoder::internal::butterfly_inv(w);
                w /= 2 * m;

                w.erase(w.begin(), w.begin() + m);
                w.resize(2 * m, 0);
                atcoder::internal::butterfly(w);
                atcoder::internal::butterfly(q);
                for (size_t i = 0; i < 2 * m; i++) {
                    w[i] *= q[i];
                }
                atcoder::internal::butterfly_inv(w);
                w /= 2 * m;

                fps df(f.derivative());
                df.resize(2 * m, 0);
                atcoder::internal::butterfly(df);
                for (size_t i = 0; i < 2 * m; i++) {
                    df[i] *= g_ft[i];
                }
                atcoder::internal::butterfly_inv(df);
                df /= 2 * m;

                std::copy(w.begin(), w.begin() + m, w.begin() + m);
                std::fill(w.begin(), w.begin() + m, 0);
                w -= df;
            }

            w = w.integral();
            for (size_t i = 0; i < 2 * m && i < n; i++) {
                w[i] += (*this)[i];
            }
            w.erase(w.begin(), w.begin() + m);
            w.resize(2 * m, 0);
            atcoder::internal::butterfly(w);
            for (size_t i = 0; i < 2 * m; i++) {
                f_ft[i] *= w[i];
            }
            atcoder::internal::butterfly_inv(f_ft);
            f_ft /= 2 * m;
            std::copy(f_ft.begin(), f_ft.begin() + m, f_ft.begin() + m);
            std::copy(f.begin(), f.begin() + m, f_ft.begin());
            f.swap(f_ft);
        }
        f.resize(n);
        return f;
    }

    fps pow(long long m) const {
        size_t n = this->size();
        size_t l = std::find_if(this->begin(), this->end(),
                                [](T x) { return x != 0; }) -
                   this->begin();
        if (l == this->size() || (l > 0 && m >= (long long)div_ceil(n, l))) {
            return fps(n, 0);
        }
        fps res(this->begin() + l, this->end());
        T c = (*this)[l];
        res /= c;
        res.resize(n, 0);
        res = (res.log() * T(m)).exp();
        res.erase(res.begin() + (n - m * l), res.end());
        res *= c.pow(m);
        std::reverse(res.begin(), res.end());
        res.resize(n, 0);
        std::reverse(res.begin(), res.end());
        return res;
    }

    void dump() {
        for (auto x : *this) printf("%d ", x.val());
        puts("");
    }

   private:
    static Factorial<T> fc;
};

template <class T>
Factorial<T> FormalPowerSeries<T>::fc;

}  // namespace rklib

#endif  // RK_FORMAL_POWER_SERIES_HPP