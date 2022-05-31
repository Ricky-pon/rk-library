#ifndef RK_AKS_PRIMALITY_TESTING_HPP
#define RK_AKS_PRIMALITY_TESTING_HPP

#include <numeric>
#include <vector>

namespace rklib {

namespace internal {

std::pair<int, bool> log2_ceil(long long n) {
    long long res = 0, tmp = 1;
    while (tmp < n) {
        tmp *= 2;
        ++res;
    }
    return {res, tmp == n};
}

long long pow_mod(long long x, long long n, long long mod) {
    long long res = 1;
    while (n) {
        if (n & 1)
            (res *= x) %= mod, --n;
        else
            (x *= x) %= mod, n >>= 1;
    }
    return res;
}

bool is_perfect_power(long long n) {
    int bmax = log2_ceil(n).first;
    for (int b = 2; b <= bmax; ++b) {
        long long low = 1, high = n;
        while (high - low > 1) {
            long long mid = (high + low) / 2;
            bool flag = false;
            long long tmp = 1;
            for (int i = 0; i < b; ++i) {
                tmp *= mid;
                if (tmp > n) {
                    flag = true;
                    break;
                }
            }
            (flag ? high : low) = mid;
        }
        if (pow_mod(low, b, n + 1) == n) return true;
    }
    return false;
}

long long sqrt_ceil(long long r) {
    for (long long i = 1;; ++i) {
        if (i * i >= r) return i;
    }
}

long long totient_function(long long r) {
    long long res = r;
    for (long long i = 2; i * i <= r; ++i) {
        if (r % i == 0) {
            res = res * (i - 1) / i;
            while (r % i == 0) r /= i;
        }
    }
    if (r > 1) res = res * (r - 1) / r;
    return res;
}

std::vector<long long> prod(std::vector<long long> &f,
                            std::vector<long long> &g, long long r,
                            long long n) {
    std::vector<long long> res(r, 0);
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < r; ++j) {
            int idx = (i + j >= r ? i + j - r : i + j);
            __int128_t tmp = __int128_t(f[i]) * __int128_t(g[j]) + res[idx];
            res[idx] = tmp % n;
        }
    }
    return res;
}

}  // namespace internal

bool is_prime(long long n) {
    // step 1
    if (internal::is_perfect_power(n)) return false;

    // step 2
    long long r = 3;
    {
        auto [threshold, eq] = internal::log2_ceil(n);
        threshold += eq;
        for (;; ++r) {
            long long tmp = 1, n_mod_r = n % r;
            int cnt = 0;
            while (true) {
                (tmp *= n_mod_r) %= r;
                ++cnt;
                if (tmp == 1 || cnt >= threshold) break;
            }
            if (cnt >= threshold) {
                break;
            }
        }
    }

    // step 3
    for (int a = 2; a < std::min(n, r + 1); ++a) {
        if (std::gcd(a, n) > 1) return false;
    }

    // step 4
    if (n <= r) return true;

    // step 5
    long long l = internal::sqrt_ceil(internal::totient_function(r)) *
                  internal::log2_ceil(n).first;
    for (int a = 1; a <= l; ++a) {
        long long m = n;
        std::vector<long long> f(r, 0), g(r, 0);
        f[0] = a % n;
        f[1] = 1;
        g[0] = 1;
        while (m) {
            if (m & 1) {
                g = internal::prod(f, g, r, n);
                --m;
            } else {
                f = internal::prod(f, f, r, n);
                m >>= 1;
            }
        }
        if (g[0] != a % n || g[n % r] != 1) return false;
        for (int i = 1; i < r; ++i)
            if (i != n % r && g[i] != 0) return false;
    }

    return true;
}

}  // namespace rklib

#endif  // RK_AKS_PRIMALITY_TESTING