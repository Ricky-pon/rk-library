#ifndef RK_SIEVE_HPP
#define RK_SIEVE_HPP

#include <algorithm>
#include <bitset>
#include <vector>

namespace rklib {

struct SimpleSieve {
    std::vector<bool> is_prime;
    std::vector<int> prime;

    SimpleSieve(int n) {
        is_prime.resize(n + 1, true);
        is_prime[0] = is_prime[1] = false;
        for (int i = 2; i <= n; ++i) {
            if (!is_prime[i]) continue;
            prime.push_back(i);
            for (int j = 2; i * j <= n; ++j) {
                is_prime[i * j] = false;
            }
        }
    }
};

struct Sieve {
    std::vector<int> min_factor, prime;

    Sieve(int n) {
        min_factor.resize(n + 1, 0);
        for (int i = 2; i <= n; ++i) {
            if (min_factor[i] == 0) {
                min_factor[i] = i;
                prime.push_back(i);
            }
            for (int x : prime) {
                if (x * i > n || x > i) break;
                min_factor[x * i] = x;
            }
        }
    }

    std::vector<std::pair<int, int>> prime_factor(int n) {
        std::vector<std::pair<int, int>> res;
        while (n > 1) {
            if (res.empty() || res.rbegin()->first != min_factor[n]) {
                res.emplace_back(min_factor[n], 1);
            } else
                ++res.rbegin()->second;
            n /= min_factor[n];
        }
        return res;
    }

    void divisor_dfs(std::vector<std::pair<int, int>> &p, int t, int cur,
                     std::vector<int> &res) {
        if (cur == (int)p.size()) {
            res.push_back(t);
            return;
        }
        divisor_dfs(p, t, cur + 1, res);
        for (int _ = 0; _ < p[cur].second; ++_) {
            t *= p[cur].first;
            divisor_dfs(p, t, cur + 1, res);
        }
    }

    std::vector<int> get_divisor(int n, bool sorted = false) {
        std::vector<int> res;
        auto p = prime_factor(n);
        divisor_dfs(p, 1, 0, res);
        if (sorted) sort(res.begin(), res.end());
        return res;
    }
};

}  // namespace rklib

#endif  // RK_SIEVE_HPP