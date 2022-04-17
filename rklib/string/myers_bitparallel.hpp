#ifndef RK_MYERS_BITPARALLEL_HPP
#define RK_MYERS_BITPARALLEL_HPP

#include <array>
#include <rklib/data_structure/bit_set.hpp>
#include <string>

namespace rklib {

template <char cmin = 'a', char cmax = 'z'>
int myers_bitparallel(std::string &s, std::string &t) {
    const size_t cnum = cmax - cmin + 1;
    size_t n = s.size(), m = t.size();

    std::array<BitSet, cnum> is;
    for (size_t c = 0; c < cnum; c++) {
        is[c] = BitSet(s.size() + 1, 0);
    }
    for (size_t i = 0; i < n; i++) {
        is[s[i] - cmin].set(i + 1, 1);
    }

    int res = n;
    BitSet r_plus(n + 1, 1), r_minus(n + 1, 0);
    r_plus.set(0, 0);
    for (size_t i = 0; i < m; i++) {
        int c = t[i] - cmin;
        BitSet d = ~((((is[c] & r_plus) + r_plus) ^ r_plus) | is[c] | r_minus);
        d.set(0, 1);
        BitSet c_plus = r_minus | ~(~d | r_plus);
        c_plus.set(0, 1);
        BitSet c_minus = ~d & r_plus;
        r_plus = (c_minus << 1) | ~(~d | (c_plus << 1));
        r_minus = ~d & (c_plus << 1);

        if (c_plus.get(n) == 1) ++res;
        if (c_minus.get(n) == 1) --res;
    }

    return res;
}

}  // namespace rklib

#endif  // RK_MYERS_BITPARALLEL_HPP