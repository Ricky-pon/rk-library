#ifndef RK_Z_ALGORITHM_HPP
#define RK_Z_ALGORITHM_HPP

#include <rklib/utility/utility.hpp>
#include <string>
#include <vector>

namespace rklib {

template <class T>
std::vector<int> z_algorithm(T& s) {
    int n = s.size();
    std::vector<int> res(n);
    res[0] = n;
    for (int i = 1, left = 0, right = 0; i < n; i++) {
        int tmp = (i + 1 <= right ? std::min(res[i - left], right - i) : 0);
        while (tmp < n && s[i + tmp] == s[tmp]) ++tmp;
        res[i] = tmp;
        if (chmax(right, i + res[i])) left = i;
    }
    return res;
}

}  // namespace rklib

#endif  // RK_Z_ALGORITHM_HPP