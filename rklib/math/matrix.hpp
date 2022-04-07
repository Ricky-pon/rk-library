#ifndef RK_MATRIX_HPP
#define RK_MATRIX_HPP

#include <assert.h>

#include <vector>

namespace rklib {

template <class T>
struct Matrix {
    std::vector<std::vector<T>> a;

    Matrix(int n, int m) : a(n, std::vector<T>(m, 0)) {}
    Matrix(int n) : a(n, std::vector<T>(n, 0)) {}

    int height() const { return a.size(); }
    int width() const { return a[0].size(); }

    inline const std::vector<T> &operator[](int i) const { return a[i]; }
    inline std::vector<T> &operator[](int i) { return a[i]; }

    static Matrix id(int n) {
        Matrix res(n);
        for (int i = 0; i < n; i++) res[i][i] = 1;
        return res;
    }

    Matrix &operator*=(const Matrix &b) {
        assert(width() == b.height());
        std::vector<std::vector<T>> c(height(), std::vector<T>(b.width()));
        for (int i = 0; i < height(); ++i) {
            for (int k = 0; k < b.height(); ++k) {
                for (int j = 0; j < b.width(); ++j) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        a.swap(c);
        return *this;
    }

    Matrix pow(long long n) const {
        auto x = (*this), res = id(height());
        while (n) {
            if (n & 1) {
                res *= x;
                --n;
            } else {
                x *= x;
                n >>= 1;
            }
        }
        return res;
    }

    Matrix operator*(const Matrix &b) const { return Matrix(*this) *= b; }

    std::vector<T> operator*(const std::vector<T> &v) {
        assert(width() == (int)v.size());
        std::vector<T> res(height(), 0);
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                res[i] += a[i][j] * v[j];
            }
        }
        return res;
    }
};

template <class T>
std::pair<int, bool> gauss_jordan(Matrix<T> &a) {
    int rnk = 0;
    bool swp = false;
    for (int j = 0; j < a.width(); ++j) {
        int pivot = -1;
        for (int i = rnk; i < a.height(); ++i) {
            if (a[i][j] != 0) {
                pivot = i;
                break;
            }
        }
        if (pivot < 0) continue;
        swap(a[pivot], a[rnk]);
        if (pivot != rnk) swp ^= true;
        for (int i = 0; i < a.height(); ++i) {
            if (i != rnk && a[i][j] != 0) {
                auto coef = a[i][j] / a[rnk][j];
                for (int k = j; k < a.width(); ++k) {
                    a[i][k] -= a[rnk][k] * coef;
                }
            }
        }
        ++rnk;
    }
    return {rnk, swp};
}

template <class T>
T determinant(Matrix<T> a) {
    auto [rnk, swp] = gauss_jordan(a);
    if (rnk < a.height()) return 0;
    T res = 1;
    for (int i = 0; i < a.height(); ++i) res *= a[i][i];
    if (swp) res = -res;
    return res;
}

template <class T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
system_of_linear_equations(Matrix<T> a, const std::vector<T> &b) {
    assert(a.height() == (int)b.size());
    Matrix<T> aug(a.height(), a.width() + 1);
    for (int i = 0; i < a.height(); ++i) {
        for (int j = 0; j < a.width(); ++j) {
            aug[i][j] = a[i][j];
        }
        aug[i][a.width()] = b[i];
    }

    auto rnk = gauss_jordan(a).first, rnk_aug = gauss_jordan(aug).first;

    std::vector<T> solution;
    std::vector<std::vector<T>> kernel;
    if (rnk < rnk_aug) return {solution, kernel};

    solution.resize(a.width(), 0);
    std::vector<bool> used(a.width(), false);
    std::vector<int> pos(rnk);
    for (int i = 0; i < rnk; ++i) {
        for (int j = 0; j < a.width(); ++j) {
            if (aug[i][j] != 0) {
                solution[j] = aug[i][a.width()] / aug[i][j];
                used[j] = true;
                pos[i] = j;
                break;
            }
        }
    }
    for (int j = 0; j < a.width(); ++j) {
        if (!used[j]) {
            std::vector<T> v(a.width(), 0);
            v[j] = 1;
            for (int i = 0; i < rnk; ++i) {
                v[pos[i]] = -aug[i][j] / aug[i][pos[i]];
            }
            kernel.push_back(v);
        }
    }
    return {solution, kernel};
}

}  // namespace rklib

#endif  // RK_MATRIX_HPP