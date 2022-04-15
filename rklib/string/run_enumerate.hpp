#ifndef RK_RUN_ENUMERATE_HPP
#define RK_RUN_ENUMERATE_HPP

#include <rklib/string/z_algorithm.hpp>
#include <vector>

namespace rklib {

namespace internal {

template <class T>
void run_enumerate_rec(int l, int r, T &s,
                       std::vector<std::tuple<int, int, int>> &runs) {
    if (r - l <= 1) return;
    int m = (l + r) / 2;
    run_enumerate_rec(l, m, s, runs);
    run_enumerate_rec(m, r, s, runs);

    auto f = [&](bool rev) {
        T t(s.begin() + l, s.begin() + r);
        if (rev) {
            std::reverse(t.begin(), t.end());
            m = l + r - m;
        }

        int len = r - l, mid = m - l;
        T tl(t.begin(), t.begin() + mid);
        std::reverse(tl.begin(), tl.end());
        T tr(t.begin() + mid, t.begin() + len);
        std::copy(t.begin(), t.end(), std::back_inserter(tr));
        auto zl = z_algorithm(tl), zr = z_algorithm(tr);
        zl.push_back(0);

        for (int k = 1; mid - k >= 0; ++k) {
            int li = m - k - zl[k], ri = m + std::min(r - m, zr[len - k]);
            if (rev) {
                std::swap(li, ri);
                li = l + r - li;
                ri = l + r - ri;
            }
            if (ri - li < 2 * k) continue;
            if (li > 0 && s[li - 1] == s[li - 1 + k]) continue;
            if (ri < int(s.size()) && s[ri] == s[ri - k]) continue;
            runs.push_back({li, ri, k});
        }
    };

    f(false);
    f(true);
}

}  // namespace internal

template <class T>
std::vector<std::tuple<int, int, int>> run_enumerate(T &s) {
    std::vector<std::tuple<int, int, int>> runs;
    internal::run_enumerate_rec(0, s.size(), s, runs);
    sort(runs.begin(), runs.end());
    std::vector<std::tuple<int, int, int>> res;
    for (int i = 0; i < (int)runs.size(); ++i) {
        if (i > 0 && std::get<0>(runs[i]) == std::get<0>(runs[i - 1]) &&
            std::get<1>(runs[i]) == std::get<1>(runs[i - 1])) {
            continue;
        }
        auto [l, r, t] = runs[i];
        res.push_back({t, l, r});
    }
    std::sort(res.begin(), res.end());
    return res;
}

}  // namespace rklib

#endif  // RK_RUN_ENUMERATE_HPP