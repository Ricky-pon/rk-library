#ifndef RK_LIFETIME_TRACKER_HPP
#define RK_LIFETIME_TRACKER_HPP

#include <cassert>
#include <map>
#include <vector>

namespace rklib {

template<class T>
struct LifetimeTracker
{
    struct Range
    {
        int l, r; // [l, r)
        T value;

        bool empty() const { return l >= r; }
    };

    struct Handle
    {
        int id;
        explicit Handle(int id)
          : id(id)
        {
        }
    };

    int q_;
    std::vector<Range> ranges_;
    std::vector<bool> open_flag_;

    explicit LifetimeTracker(int q)
      : q_(q)
    {
    }

    Handle open(int t, const T& value)
    {
        Handle h{ int(ranges_.size()) };
        ranges_.push_back(Range(t, q_, value));
        open_flag_.push_back(true);
        return h;
    }

    void close(int t, Handle h)
    {
        assert(open_flag_[h.id]);
        ranges_[h.id].r = t;
        open_flag_[h.id] = false;
    }

    std::vector<Range> build()
    {
        std::vector<Range> res;
        for (auto rg : ranges_) {
            if (!rg.empty()) {
                res.push_back(rg);
            }
        }
        return res;
    }
};

template<class T, class Key = T>
struct KeyedLifetimeTracker
{
    LifetimeTracker<T> tracker_;
    std::map<Key, typename LifetimeTracker<T>::Handle> handles_;

    explicit KeyedLifetimeTracker(int q)
      : tracker_(q)
    {
    }

    void open(int t, const Key& key, const T& value)
    {
        handles_.insert_or_assign(key, tracker_.open(t, value));
    }

    void close(int t, const Key& key)
    {
        tracker_.close(t, handles_.at(key));
        handles_.erase(key);
    }

    std::vector<typename LifetimeTracker<T>::Range> build()
    {
        return tracker_.build();
    }
};

}

#endif // RK_LIFETIME_TRACKER_HPP
