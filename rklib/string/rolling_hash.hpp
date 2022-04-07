#ifndef RK_ROLLING_HASH_HPP
#define RK_ROLLING_HASH_HPP

#include <string>
#include <vector>

namespace rklib {

struct RollingHash {
   private:
    using ulint = unsigned long long;
    const ulint hash_mod = (1ULL << 61ULL) - 1, hash_base = 1'000'000'007;
    const ulint mask_30bit = (1ULL << 30ULL) - 1,
                mask_31bit = (mask_30bit << 1ULL) + 1;
    const ulint num = hash_mod << 2LL;
    std::vector<ulint> pow_table, hash;
    int n;

    ulint calc_mul31(ulint a, ulint b) {
        ulint au = a >> 31ULL, ad = a & mask_31bit;
        ulint mid = au * b, midu = mid >> 30ULL, midd = mid & mask_30bit;
        return midu + (midd << 31ULL) + ad * b;
    }

    ulint calc_mul61(ulint a, ulint b) {
        ulint au = a >> 31ULL, ad = a & mask_31bit;
        ulint bu = b >> 31ULL, bd = b & mask_31bit;
        ulint mid = au * bd + ad * bu, midu = mid >> 30ULL,
              midd = mid & mask_30bit;
        return 2 * au * bu + midu + (midd << 31ULL) + ad * bd;
    }

    ulint calc_mod(ulint x) {
        x = (x & hash_mod) + (x >> 61ULL);
        if (x > hash_mod) x -= hash_mod;
        return x;
    }

   public:
    RollingHash(std::string &s) {
        n = s.size();
        pow_table.resize(n + 1);
        hash.resize(n + 1);
        pow_table[0] = 1ULL;
        for (int i = 0; i < n; ++i)
            pow_table[i + 1] = calc_mod(calc_mul61(pow_table[i], hash_base));
        hash[0] = 0ULL;
        for (int i = 0; i < n; ++i)
            hash[i + 1] = calc_mod(calc_mul61(hash[i], hash_base) + s[i]);
    }

    ulint slice(int l, int r) {
        return calc_mod(hash[r] + num - calc_mul61(hash[l], pow_table[r - l]));
    }
};

}  // namespace rklib

#endif  // RK_ROLLING_HASH_HPP