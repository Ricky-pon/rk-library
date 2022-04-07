#ifndef RK_DISJOINT_SPARSE_TABLE_HPP
#define RK_DISJOINT_SPARSE_TABLE_HPP

#include <rklib/utility/utility.hpp>
#include <vector>

namespace rklib {

template <class S, S (*op)(S, S), S (*e)()>
struct DisjointSparseTable {
   public:
    DisjointSparseTable() : DisjointSparseTable(1) {}
    DisjointSparseTable(int n) : DisjointSparseTable(std::vector<S>(n, e())) {}
    DisjointSparseTable(std::vector<S> &v) : w((int)v.size()) {
        if (w == 0) {
            h = 0;
            return;
        }
        h = (w == 1 ? 1 : 32 - __builtin_clz(w - 1));
        table.resize(h);
        for (int i = 0; i < h; i++) {
            table[i].resize(w);

            int step = 1 << (i + 1);
            for (int l = 0, s = (1 << i) - 1, t = 1 << i,
                     r = (1 << (i + 1)) - 1;
                 l < w; l += step, s += step, t += step, r += step) {
                chmin(s, w - 1);
                table[i][s] = v[s];
                for (int j = s - 1; j >= l; j--) {
                    table[i][j] = op(v[j], table[i][j + 1]);
                }
                if (s == w - 1) break;

                chmin(r, w - 1);
                table[i][t] = v[t];
                for (int j = t + 1; j < r + 1; j++) {
                    table[i][j] = op(table[i][j - 1], v[j]);
                }
                if (r == w - 1) break;
            }
        }
    }

    S prod(int l, int r) {
        if (l == r) return e();
        --r;
        if (l == r) return table[0][l];
        int t = 31 - __builtin_clz(l ^ r);
        return op(table[t][l], table[t][r]);
    }

   private:
    int h, w;
    std::vector<std::vector<S>> table;
};

}  // namespace rklib

#endif  // RK_DISJOINT_SPARSE_TABLE_HPP