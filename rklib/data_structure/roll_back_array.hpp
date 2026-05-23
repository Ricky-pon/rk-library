#ifndef RK_ROLL_BACK_ARRAY_HPP
#define RK_ROLL_BACK_ARRAY_HPP

#include <cassert>
#include <vector>

namespace rklib {

template<class T>
struct RollbackArray
{
    int n_;
    std::vector<T> data_;
    std::vector<std::pair<int, T>> history_; // {idx, old_val}

    // === 構築 ===

    RollbackArray(std::vector<T> v)
      : n_(v.size())
      , data_(v)
    {
    }

    // === 配列操作 ===

    T get(int i)
    {
        assert(0 <= i && i < n_);
        return data_[i];
    }

    void set(int i, T val)
    {
        assert(0 <= i && i < n_);
        history_.emplace_back(i, data_[i]);
        data_[i] = val;
    }

    // === 履歴操作 ===

    struct RollbackScope
    {
        RollbackArray& ra_;
        int stamp_;

        RollbackScope(RollbackArray& ra)
          : ra_(ra)
          , stamp_(int(ra.history_.size()))
        {
        }

        ~RollbackScope()
        {
            while (int(ra_.history_.size()) > stamp_) {
                auto [i, v] = ra_.history_.back();
                ra_.history_.pop_back();
                ra_.data_[i] = v;
            }
        }

        RollbackScope(RollbackScope&&) = delete;
        RollbackScope& operator=(RollbackScope&&) = delete;
    };

    RollbackScope rollback_scope() { return RollbackScope(*this); }
};

}

#endif // RK_ROLL_BACK_ARRAY_HPP
