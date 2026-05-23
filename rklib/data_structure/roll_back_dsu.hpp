#ifndef RK_ROLL_BACK_DSU_HPP
#define RK_ROLL_BACK_DSU_HPP

#include "rklib/data_structure/roll_back_array.hpp"

namespace rklib {

struct RollbackDSU
{
    RollbackArray<int> data_; // parent or -size

    RollbackDSU(int n)
      : data_(std::vector<int>(n, -1))
    {
    }

    int leader(int v)
    {
        while (data_.get(v) >= 0) {
            v = data_.get(v);
        }
        return v;
    }

    int size(int v) { return -data_.get(leader(v)); }

    bool same(int x, int y) { return leader(x) == leader(y); }

    bool merge(int x, int y)
    {
        x = leader(x);
        y = leader(y);
        if (x == y) {
            return false;
        }
        if (size(x) < size(y)) {
            std::swap(x, y);
        }
        data_.set(x, data_.get(x) + data_.get(y));
        data_.set(y, x);
        return true;
    }

    using RollbackScope = RollbackArray<int>::RollbackScope;
    RollbackScope rollback_scope() { return data_.rollback_scope(); }
};

}

#endif // RK_ROLL_BACK_DSU_HPP
