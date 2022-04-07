#ifndef RK_FACTORIAL_HPP
#define RK_FACTORIAL_HPP

#include <atcoder/modint>
#include <vector>

namespace rklib {

template <class T>
struct Factorial {
   public:
    Factorial() : Factorial(0) {}
    Factorial(int n) {
        fc.resize(n + 1);
        inv_fc.resize(n + 1);
        fc[0] = 1;
        for (int i = 0; i < n; ++i) fc[i + 1] = fc[i] * (i + 1);

        inv_fc[n] = 1 / fc[n];
        for (int i = n - 1; i >= 0; --i) inv_fc[i] = inv_fc[i + 1] * (i + 1);
    }

    T fact(int n) {
        if (n >= (int)fc.size()) extend(n);
        return fc[n];
    }

    T inv_fact(int n) {
        if (n >= (int)fc.size()) extend(n);
        return inv_fc[n];
    }

    T inv(int n) {
        assert(n > 0);
        if (n >= (int)fc.size()) extend(n);
        return inv_fc[n] * fc[n - 1];
    }

    T comb(int n, int r) {
        if (n < r || r < 0) return 0;
        if (n >= (int)fc.size()) extend(n);
        return fc[n] * inv_fc[r] * inv_fc[n - r];
    }

    T perm(int n, int r) {
        if (n < r || r < 0) return 0;
        if (n >= (int)fc.size()) extend(n);
        return fc[n] * inv_fc[n - r];
    }

   private:
    std::vector<T> fc;
    std::vector<T> inv_fc;

    void extend(int n) {
        int l = fc.size();
        int r = l;
        while (r <= n) r *= 2;

        fc.resize(r);
        inv_fc.resize(r);

        for (int i = l; i < r; ++i) fc[i] = fc[i - 1] * i;

        inv_fc[r - 1] = 1 / fc[r - 1];
        for (int i = r - 2; i >= l; --i) inv_fc[i] = inv_fc[i + 1] * (i + 1);
    }
};

}  // namespace rklib

#endif  // RK_FACTORIAL_HPP