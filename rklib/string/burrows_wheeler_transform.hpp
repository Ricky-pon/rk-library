#ifndef RK_BURROWS_WHEELER_TRANSFORM_HPP
#define RK_BURROWS_WHEELER_TRANSFORM_HPP

#include <atcoder/string>
#include <rklib/data_structure/bit_vector.hpp>
#include <string>
#include <vector>

namespace rklib {

template <char cmin = 'a', char cmax = 'z'>
struct BurrowsWheelerTransform {
   public:
    BurrowsWheelerTransform(std::string &s, size_t step)
        : n(s.size() + 1), step(step) {
        std::vector<int> v(s.size());
        for (size_t i = 0; i < s.size(); i++) {
            v[i] = (s[i] - cmin) + 1;
        }
        sa = atcoder::suffix_array(v);
        v.push_back(0);
        sa.insert(sa.begin(), n - 1);
        bwt.resize(n);
        for (size_t i = 0; i < n; ++i) {
            bwt[i] = (sa[i] == 0 ? v[n - 1] : v[sa[i] - 1]);
        }

        cnt_smaller.resize(cnum, 0);
        cnt.resize(cnum);
        for (size_t c = 0; c < cnum; ++c) cnt[c].resize((n + 1) / step + 1, 0);
        std::vector<int> table(cnum, 0);
        for (size_t i = 0; i < n; ++i) {
            if (bwt[i] + 1 < (int)cnum) ++cnt_smaller[bwt[i] + 1];
            ++table[bwt[i]];
            if ((i + 1) % step == 0) {
                for (size_t c = 0; c < cnum; ++c)
                    cnt[c][(i + 1) / step] = table[c];
            }
        }
        std::partial_sum(cnt_smaller.begin(), cnt_smaller.end(),
                         cnt_smaller.begin());
    }

    std::pair<int, int> fm_index(std::string &s) {
        int l = 0, r = n;
        for (int i = (int)s.size() - 1; i >= 0; --i) {
            int c = (s[i] - cmin) + 1;
            l = cnt_smaller[c] + get_cnt(c, l);
            r = cnt_smaller[c] + get_cnt(c, r);
        }
        return {l, r};
    }

    std::vector<int> find_all(std::string &s) {
        auto [l, r] = fm_index(s);
        std::vector<int> res(sa.begin() + l, sa.begin() + r);
        std::sort(res.begin(), res.end());
        return res;
    }

   private:
    const size_t cnum = (cmax - cmin) + 2;
    size_t n, step;
    std::vector<int> bwt, cnt_smaller, sa;
    std::vector<std::vector<int>> cnt;

    int get_cnt(int c, int k) {
        int ret = cnt[c][k / step];
        for (int i = k / step * step; i < k; ++i) ret += (bwt[i] == c);
        return ret;
    }
};

template <char cmin = 'a', char cmax = 'z'>
struct BurrowsWheelerTransformBitVector {
   public:
    BurrowsWheelerTransformBitVector(std::string &s, size_t step)
        : n(s.size() + 1), step(step) {
        std::vector<int> v(s.size());
        for (size_t i = 0; i < s.size(); i++) {
            v[i] = (s[i] - cmin) + 1;
        }
        sa = atcoder::suffix_array(v);
        v.push_back(0);
        sa.insert(sa.begin(), n - 1);
        bwt.resize(n);
        for (size_t i = 0; i < n; ++i) {
            bwt[i] = (sa[i] == 0 ? v[n - 1] : v[sa[i] - 1]);
        }

        vs.resize(cnum);
        for (size_t c = 1; c < cnum; c++) {
            std::vector<int> b(n, 0);
            for (size_t i = 0; i < n; i++) {
                if (bwt[i] == int(c)) b[i] = 1;
            }
            vs[c] = {b};
        }

        cnt_smaller.resize(cnum, 0);
        for (size_t i = 0; i < n; ++i) {
            if (bwt[i] + 1 < (int)cnum) ++cnt_smaller[bwt[i] + 1];
        }
        std::partial_sum(cnt_smaller.begin(), cnt_smaller.end(),
                         cnt_smaller.begin());
    }

    std::pair<int, int> fm_index(std::string &s) {
        int l = 0, r = n;
        for (int i = (int)s.size() - 1; i >= 0; --i) {
            int c = (s[i] - cmin) + 1;
            l = cnt_smaller[c] + vs[c].rank(l, 1);
            r = cnt_smaller[c] + vs[c].rank(r, 1);
        }
        return {l, r};
    }

    std::vector<int> find_all(std::string &s) {
        auto [l, r] = fm_index(s);
        std::vector<int> res(sa.begin() + l, sa.begin() + r);
        std::sort(res.begin(), res.end());
        return res;
    }

   private:
    const size_t cnum = (cmax - cmin) + 2;
    size_t n, step;
    std::vector<int> bwt, cnt_smaller, sa;
    std::vector<BitVector<int>> vs;
};

}  // namespace rklib

#endif  // RK_BURROWS_WHEELER_TRANSFORM_HPP