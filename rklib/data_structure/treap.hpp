#ifndef RK_TREAP_HPP
#define RK_TREAP_HPP

#include <vector>

namespace rklib {

template <class S, S (*op)(S, S), S (*e)(), class F, S (*mapping)(F, S),
          F (*composition)(F, F), F (*id)()>
struct Treap {
   public:
    Treap() : Treap(0) {}
    Treap(int n) : Treap(std::vector<S>(n, e())) {}
    Treap(std::vector<S> &v) : root(nullptr) { root = build(0, v.size(), v); }

    struct Node {
        Node *l, *r;
        int pri, cnt;
        S val, sum;
        F lazy;
        bool rev;

        Node(S val = e())
            : l(nullptr),
              r(nullptr),
              cnt(1),
              val(val),
              sum(val),
              lazy(id()),
              rev(false) {}
    };

    void insert(int p, S x) {
        auto [l, r] = split(root, p);
        auto m = new Node(x);
        m->pri = gen();
        root = merge(merge(l, m), r);
    }

    void erase(int p) {
        auto [lm, r] = split(root, p + 1);
        auto [l, m] = split(lm, p);
        free(m);
        root = merge(l, r);
    }

    void reverse(int l, int r) {
        auto [lm_node, r_node] = split(root, r);
        auto [l_node, m_node] = split(lm_node, l);
        std::swap(m_node->l, m_node->r);
        m_node->rev ^= true;
        push(m_node);
        root = merge(merge(l_node, m_node), r_node);
    }

    S prod(int l, int r) {
        auto [lm_node, r_node] = split(root, r);
        auto [l_node, m_node] = split(lm_node, l);
        auto ret = m_node->sum;
        root = merge(merge(l_node, m_node), r_node);
        return ret;
    }

    void apply(int l, int r, F f) {
        auto [lm_node, r_node] = split(root, r);
        auto [l_node, m_node] = split(lm_node, l);
        all_apply(m_node, f);
        root = merge(merge(l_node, m_node), r_node);
    }

    void shift(int l, int r) {
        auto [lmn_node, r_node] = split(root, r);
        auto [lm_node, n_node] = split(lmn_node, r - 1);
        auto [l_node, m_node] = split(lm_node, l);
        root = merge(merge(merge(l_node, n_node), m_node), r_node);
    }

   private:
    Node *root;

    inline int gen() {
        static int x = 123456789;
        static int y = 362436069;
        static int z = 521288629;
        static int w = 88675123;
        int t;

        t = x ^ (x << 11);
        x = y;
        y = z;
        z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }

    int count(Node *t) { return t ? t->cnt : 0; }

    S sum(Node *t) { return t ? t->sum : e(); }

    Node *update(Node *t) {
        t->cnt = count(t->l) + count(t->r) + 1;
        t->sum = op(op(sum(t->l), t->val), sum(t->r));
        return t;
    }

    void all_apply(Node *t, F f) {
        if (!t) return;
        t->val = mapping(f, t->val);
        t->sum = mapping(f, t->sum);
        t->lazy = composition(f, t->lazy);
    }

    void toggle(Node *t) {
        if (!t) return;
        std::swap(t->l, t->r);
        t->rev ^= true;
        // swap(sum, sum_rev);
    }

    void push(Node *t) {
        all_apply(t->l, t->lazy);
        all_apply(t->r, t->lazy);
        t->lazy = id();

        if (t->rev) {
            toggle(t->l);
            toggle(t->r);
            t->rev = false;
        }
        update(t);
    }

    Node *merge(Node *l, Node *r) {
        if (!l || !r) return l ? l : r;
        if (l->pri > r->pri) {
            push(l);
            l->r = merge(l->r, r);
            return update(l);
        } else {
            push(r);
            r->l = merge(l, r->l);
            return update(r);
        }
    }

    std::pair<Node *, Node *> split(Node *t, int k) {
        if (!t) return {t, t};
        push(t);
        if (k <= count(t->l)) {
            auto [l, r] = split(t->l, k);
            t->l = r;
            return {l, update(t)};
        } else {
            auto [l, r] = split(t->r, k - count(t->l) - 1);
            t->r = l;
            return {update(t), r};
        }
    }

    Node *build(int l, int r, std::vector<S> &v) {
        if (l == (int)v.size()) return nullptr;
        if (r - l == 1) {
            auto t = new Node(v[l]);
            t->pri = gen();
            return t;
        }
        return merge(build(l, (l + r) / 2, v), build((l + r) / 2, r, v));
    }
};

}  // namespace rklib

#endif  // RK_TREAP_HPP