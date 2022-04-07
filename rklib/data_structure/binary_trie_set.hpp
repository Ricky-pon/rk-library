#ifndef RK_BINARY_TRIE_SET_HPP
#define RK_BINARY_TRIE_SET_HPP

#include <utility>

namespace rklib {

template <typename T, int len = 31>
struct BinaryTrieSet {
   public:
    BinaryTrieSet() : root(new Node()) {}

    int size() { return count_child(root); }

    void insert(T val) { _insert(root, val, len - 1); }

    void erase(T val) { _erase(root, val, len - 1); }

    int count(T x) {
        auto t = root;
        for (int i = len - 1; i >= 0; --i) {
            if (!t) return 0;
            push(t, i);
            if (!(x >> i & 1))
                t = t->left;
            else
                t = t->right;
        }
        return count_child(t);
    }

    T get_kth(int k) {
        ++k;
        T ret = 0;
        auto t = root;
        for (int i = len - 1; i >= 0; --i) {
            push(t, i);
            if (count_child(t->left) >= k) {
                t = t->left;
                ret <<= 1;
            } else {
                k -= count_child(t->left);
                t = t->right;
                ret <<= 1;
                ++ret;
            }
        }
        return ret;
    }

    void apply_xor(T x) { root->lz_xor ^= x; }

   private:
    struct Node {
        Node *left, *right;
        int cnt;
        T lz_xor;
        Node() : left(nullptr), right(nullptr), cnt(0), lz_xor(0) {}
    };

    Node *root;

    int count_child(Node *t) { return t ? t->cnt : 0; }

    void push(Node *t, int dep) {
        if (t->lz_xor >> dep & 1) std::swap(t->left, t->right);
        if (t->left) t->left->lz_xor ^= t->lz_xor;
        if (t->right) t->right->lz_xor ^= t->lz_xor;
        t->lz_xor = 0;
    }

    Node *_insert(Node *t, T val, int dep) {
        if (!t) t = new Node();
        if (dep == -1) {
            t->cnt = 1;
            return t;
        }
        push(t, dep);
        if (!(val >> dep & 1))
            t->left = _insert(t->left, val, dep - 1);
        else
            t->right = _insert(t->right, val, dep - 1);
        t->cnt = count_child(t->left) + count_child(t->right);
        return t;
    }

    Node *_erase(Node *t, T val, int dep) {
        if (count_child(t) == 0) return t;
        if (dep == -1) {
            t->cnt = 0;
            return t;
        }
        push(t, dep);
        if (!(val >> dep & 1))
            t->left = _erase(t->left, val, dep - 1);
        else
            t->right = _erase(t->right, val, dep - 1);
        t->cnt = count_child(t->left) + count_child(t->right);
        return t;
    }
};

}  // namespace rklib

#endif  // RK_BINARY_TRIE_SET_HPP