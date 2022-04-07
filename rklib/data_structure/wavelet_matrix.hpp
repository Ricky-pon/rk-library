#ifndef RK_WAVELET_MATRIX_HPP
#define RK_WAVELET_MATRIX_HPP

#include <array>
#include <rklib/utility/utility.hpp>
#include <vector>

namespace rklib {

template <class T>
struct BitVector {
   public:
    BitVector() : BitVector(0) {}
    BitVector(int n) : BitVector(std::vector<T>(n, 0), 0) {}
    BitVector(const std::vector<T> &a, int d) {
        int num = (a.size() >> lg) + 1;
        sum.resize(num + 1, 0);
        bit.resize(num, 0);
        for (int i = 0; i < (int)a.size(); i++) {
            if (a[i] >> d & 1) bit[i >> lg] |= (uint64_t)1 << (i & (w - 1));
        }
        for (int k = 0; k < num; k++) {
            sum[k + 1] = sum[k] + __builtin_popcountll(bit[k]);
        }
    }

    int rank(int i, int x) {
        int k = i >> lg, l = i & (w - 1);
        int res =
            sum[k] + __builtin_popcountll(bit[k] & (((uint64_t)1 << l) - 1));
        return x == 1 ? res : i - res;
    }

   private:
    static constexpr int lg = 6;
    static constexpr int w = 1 << lg;
    std::vector<int> sum;
    std::vector<uint64_t> bit;
};

template <class T, int w = 31>
struct WaveletMatrix {
   public:
    WaveletMatrix() : WaveletMatrix(0) {}
    WaveletMatrix(int n) : WaveletMatrix(std::vector<T>(n, 0)) {}
    WaveletMatrix(std::vector<T> a) : n(a.size()) {
        std::vector<std::vector<T>> b(2, std::vector<T>(n));
        copy(a.begin(), a.end(), b[(w + 1) & 1].begin());
        for (int i = w - 1; i >= 0; i--) {
            v[i] = {b[i & 1], i};
            int l = 0, r = n - 1;
            for (int j = 0; j < n; j++) {
                if ((b[i & 1][j] >> i) & 1)
                    b[(i + 1) & 1][r--] = b[i & 1][j];
                else
                    b[(i + 1) & 1][l++] = b[i & 1][j];
            }
            reverse(b[(i + 1) & 1].begin() + l, b[(i + 1) & 1].end());
        }
    }

    T quantile(int l, int r, int k) {
        T res = 0;
        for (int i = w - 1; i >= 0; i--) {
            int lz = v[i].rank(l, 0), rz = v[i].rank(r, 0), cntz = rz - lz;
            if (cntz > k) {
                l = lz;
                r = rz;
            } else {
                res |= T(1) << i;
                k -= cntz;
                int t = v[i].rank(n, 0);
                l = t + l - lz;
                r = t + r - rz;
            }
        }
        return res;
    }

    int freq(int l, int r, T lower, T upper) {
        return _freq(l, r, upper) - _freq(l, r, lower);
    }

    int freq(int l, int r, T x) {
        for (int i = w - 1; i >= 0; i--) {
            int lz = v[i].rank(l, 0), rz = v[i].rank(r, 0);
            if (!(x >> i & 1)) {
                l = lz;
                r = rz;
            } else {
                int t = v[i].rank(n, 0);
                l = t + l - lz;
                r = t + r - rz;
            }
        }
        return r - l;
    }

   private:
    int n;
    std::array<BitVector<T>, w> v;

    int _freq(int l, int r, T upper) {
        int res = 0;
        for (int i = w - 1; i >= 0; i--) {
            int lz = v[i].rank(l, 0), rz = v[i].rank(r, 0), cntz = rz - lz;
            if (!(upper >> i & 1)) {
                l = lz;
                r = rz;
            } else {
                res += cntz;
                int t = v[i].rank(n, 0);
                l = t + l - lz;
                r = t + r - rz;
            }
        }
        return res;
    }
};

}  // namespace rklib

#endif  // RK_WAVELET_MATRIX_HPP