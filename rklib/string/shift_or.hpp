#ifndef RK_SHIFT_OR_HPP
#define RK_SHIFT_OR_HPP

#include <array>
#include <rklib/data_structure/bit_set.hpp>
#include <string>

namespace rklib {

template <char cmin = 'a', char cmax = 'z'>
struct ShiftOr {
    const int cnum = cmax - cmin + 1;
    size_t n;
    std::vector<BitSet> bs;

    ShiftOr(std::string &s) : n(s.size()) {
        bs.resize(cnum, BitSet(n, 1));
        for (size_t i = 0; i < n; i++) {
            bs[s[i] - cmin].set(i, 0);
        }
    }

    std::vector<int> find(std::string &s) {
        std::vector<int> res;
        BitSet v(n, 1);
        for (size_t i = 0; i < s.size(); i++) {
            v <<= 1;
            v |= bs[s[i] - cmin];
            if (v.get(n - 1) == 0) {
                res.push_back(i + 1 - n);
            }
        }
        return res;
    }
};

}  // namespace rklib

#endif  // RK_SHIFT_OR_HPP