#ifndef RK_BITVECTOR_HPP
#define RK_BITVECTOR_HPP

#include <cstdint>
#include <vector>

namespace rklib {

template <class T>
struct BitVector {
   public:
    BitVector() : BitVector(0) {}
    BitVector(int n) : BitVector(std::vector<T>(n, 0), 0) {}
    BitVector(const std::vector<T> &a, int d = 0) {
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

}  // namespace rklib

#endif  // RK_BITVECTOR_HPP