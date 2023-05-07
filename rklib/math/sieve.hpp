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
   public:
    Sieve(int n) {
        mpf.resize(n + 1, 0);
        mpf[1] = 1;
        for (int i = 2; i <= n; ++i) {
            if (mpf[i] == 0) {
                mpf[i] = i;
                ps.push_back(i);
            }
            for (int p : ps) {
                if (p * i > n || p > mpf[i]) break;
                mpf[p * i] = p;
            }
        }
    }

    bool is_prime(int n) { return mpf[n] == n; }

    int min_prime_factor(int n) { return mpf[n]; }

    std::vector<std::pair<int, int>> prime_factorize(int n) {
        std::vector<std::pair<int, int>> res;
        while (n > 1) {
            if (res.empty() || res.rbegin()->first != mpf[n]) {
                res.emplace_back(mpf[n], 1);
            } else
                ++res.rbegin()->second;
            n /= mpf[n];
        }
        return res;
    }

    std::vector<int> divisors(int n, bool sorted = false) {
        std::vector<int> res;
        auto p = prime_factorize(n);
        divisor_dfs(p, 1, 0, res);
        if (sorted) sort(res.begin(), res.end());
        return res;
    }

    std::vector<int> primes() { return ps; }

   private:
    std::vector<int> mpf, ps;

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
};

}  // namespace rklib

#endif  // RK_SIEVE_HPP