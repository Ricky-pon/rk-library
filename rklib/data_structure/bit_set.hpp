#ifndef RK_BIT_SET_HPP
#define RK_BIT_SET_HPP

#include <rklib/utility/utility.hpp>
#include <vector>

namespace rklib {

struct BitSet {
    using ulint = unsigned long long;

    size_t n;
    std::vector<ulint> bit;
    static constexpr size_t lg = 6, w = 64, mod_w = w - 1;
    static constexpr ulint mask = (((1ULL << (w - 1)) - 1) << 1) + 1;

    BitSet() : BitSet(0) {}
    BitSet(size_t n) : BitSet(n, 0) {}
    BitSet(size_t n, int v) : n(n) {
        ulint val = (v == 0 ? 0ULL : mask);
        bit.resize(div_ceil(n, w), val);
        size_t r = n & mod_w;
        if (r > 0) *bit.rbegin() &= mask >> (w - r);
    }
    BitSet(const std::vector<int>& v) : n(v.size()) {
        bit.resize(div_ceil(n, w), 0ULL);
        for (size_t i = 0; i < bit.size(); i++) {
            for (size_t j = 0; j < w && i * w + j < n; j++) {
                if (v[i * w + j] != 0) bit[i] |= 1ULL << j;
            }
        }
    }

    void set(size_t pos, int val) {
        assert(0 <= pos && pos < n);
        size_t i = pos >> lg, j = pos & mod_w;
        ulint left = (j > 0 ? mask >> (w - j) : 0ULL);
        ulint right = (j + 1 < w ? (mask >> (j + 1)) << (j + 1) : 0);
        bit[i] = (bit[i] & left) | (bit[i] & right);
        if (val != 0) bit[i] |= 1ULL << j;
    }

    int get(size_t pos) {
        assert(0 <= pos && pos < n);
        size_t i = pos >> lg, j = pos & mod_w;
        return bit[i] >> j & 1;
    }

    size_t size() { return n; }

    BitSet& operator|=(const BitSet& rhs) {
        if (this->bit.size() < rhs.bit.size())
            this->bit.resize(rhs.bit.size(), 0);
        for (size_t i = 0; i < this->bit.size() && i < rhs.bit.size(); i++) {
            this->bit[i] |= rhs.bit[i];
        }
        return *this;
    }
    BitSet& operator&=(const BitSet& rhs) {
        if (this->bit.size() < rhs.bit.size())
            this->bit.resize(rhs.bit.size(), 0);
        for (size_t i = 0; i < this->bit.size() && i < rhs.bit.size(); i++) {
            this->bit[i] &= rhs.bit[i];
        }
        return *this;
    }
    BitSet& operator^=(const BitSet& rhs) {
        if (this->bit.size() < rhs.bit.size())
            this->bit.resize(rhs.bit.size(), 0);
        for (size_t i = 0; i < this->bit.size() && i < rhs.bit.size(); i++) {
            this->bit[i] ^= rhs.bit[i];
        }
        return *this;
    }
    BitSet& operator+=(const BitSet& rhs) {
        if (this->bit.size() < rhs.bit.size())
            this->bit.resize(rhs.bit.size(), 0);
        ulint carry = 0ULL;
        for (size_t i = 0; i < this->bit.size() && i < rhs.bit.size(); i++) {
            ulint tmp = 0;
            bool flag =
                __builtin_uaddll_overflow(this->bit[i], rhs.bit[i], &tmp);
            flag |= __builtin_uaddll_overflow(tmp, carry, &tmp);
            this->bit[i] = tmp;
            carry = (flag ? 1ULL : 0ULL);
        }
        return *this;
    }
    BitSet& operator<<=(const ulint& rhs) {
        if (rhs >= n) {
            std::fill(this->bit.begin(), this->bit.end(), 0ULL);
            return *this;
        }
        size_t q = rhs >> lg, r = rhs & mod_w;
        if (q > 0) {
            for (int i = int(bit.size()) - 1; i >= 0; i--) {
                this->bit[i] = (i >= int(q) ? this->bit[i - q] : 0ULL);
            }
        }
        if (r > 0) {
            for (int i = int(bit.size()) - 1; i >= 1; i--) {
                this->bit[i] =
                    (this->bit[i - 1] >> (w - r)) | (this->bit[i] << r);
            }
            this->bit[0] <<= r;
        }
        return *this;
    }

    BitSet operator~() const {
        BitSet res(n);
        for (size_t i = 0; i < this->bit.size(); i++) {
            res.bit[i] = ~this->bit[i];
        }
        size_t j = n & mod_w;
        if (j > 0) {
            *res.bit.rbegin() &= mask >> (w - j);
        }
        return res;
    }

    friend BitSet operator|(const BitSet& lhs, const BitSet& rhs) {
        return BitSet(lhs) |= rhs;
    }
    friend BitSet operator&(const BitSet& lhs, const BitSet& rhs) {
        return BitSet(lhs) &= rhs;
    }
    friend BitSet operator^(const BitSet& lhs, const BitSet& rhs) {
        return BitSet(lhs) ^= rhs;
    }
    friend BitSet operator+(const BitSet& lhs, const BitSet& rhs) {
        return BitSet(lhs) += rhs;
    }
    friend BitSet operator<<(const BitSet& lhs, const ulint& rhs) {
        return BitSet(lhs) <<= rhs;
    }

   private:
};

}  // namespace rklib

#endif  // RK_BIT_SET_HPP