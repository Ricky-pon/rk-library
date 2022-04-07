#ifndef RK_DUAL_SEGTREE_HPP
#define RK_DUAL_SEGTREE_HPP

#include <cassert>
#include <vector>

namespace rklib {

template <class S, S (*op)(S, S), S (*e)()>
struct DualSegTree {
   public:
    DualSegTree() : DualSegTree(0) {}
    DualSegTree(int n) : DualSegTree(std::vector<S>(n, e())) {}
    DualSegTree(const std::vector<S> &v) : _n(int(v.size())) {
        size = 1;
        log = 0;
        while (size < _n) size <<= 1, ++log;
        d = std::vector<S>(2 * size, e());
        for (int i = 0; i < _n; i++) d[size + i] = v[i];
    }

    S get(int p) {
        assert(0 <= p && p < _n);
        p += size;
        S sum = e();
        for (int i = 0; i <= log; i++) sum = op(sum, d[p >> i]);
        return sum;
    }

    void prod(int l, int r, S x) {
        assert(0 <= l && l <= r && r <= _n);
        l += size;
        r += size;

        while (l < r) {
            if (l & 1) d[l] = op(d[l], x), ++l;
            if (r & 1) --r, d[r] = op(d[r], x);
            l >>= 1;
            r >>= 1;
        }
    }

    void all_prod(S x) { d[1] = op(d[1], x); }

   private:
    int _n, size, log;
    std::vector<S> d;
};

}  // namespace rklib

#endif  // RK_DUAL_SEGTREE_HPP