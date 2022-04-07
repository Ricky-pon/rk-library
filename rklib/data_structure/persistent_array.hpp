#ifndef RK_PERSISTENT_ARRAY_HPP
#define RK_PERSISTENT_ARRAY_HPP

#include <unordered_map>
#include <vector>

namespace rklib {

template <typename T>
struct PersistentArray {
    struct Node {
        T val;
        Node *ch[2];

        Node() : val{}, ch{nullptr, nullptr} {}
    };

    std::unordered_map<int, Node *> roots;

    PersistentArray() {}

    PersistentArray(std::vector<T> &a) { build(a); }

    PersistentArray(int n, T val) {
        std::vector<T> a(n, val);
        build(a);
    }

    void build(std::vector<T> &a) {
        Node *root = new Node();
        roots[-1] = root;
        for (int i = 0; i < (int)a.size(); i++) {
            Node *node = root;
            int idx = i;
            while (true) {
                if (idx == 0) {
                    node->val = a[i];
                    break;
                }
                if (!node->ch[idx & 1]) node->ch[idx & 1] = new Node();
                node = node->ch[idx & 1];
                idx >>= 1;
            }
        }
    }

    void set(int now, int prv, int idx, T val) {
        if (roots.find(now) == roots.end()) {
            roots[now] = new Node();
            roots[now]->val = roots[prv]->val;
            for (int c = 0; c < 2; c++) roots[now]->ch[c] = roots[prv]->ch[c];
        }
        Node *now_node = roots[now], *prv_node = roots[prv];
        while (true) {
            if (idx == 0) {
                now_node->val = val;
                return;
            }
            if (!now_node->ch[idx & 1] ||
                now_node->ch[idx & 1] == prv_node->ch[idx & 1]) {
                now_node->ch[idx & 1] = new Node();
                now_node->ch[idx & 1]->val = prv_node->ch[idx & 1]->val;
                for (int c = 0; c < 2; c++)
                    now_node->ch[idx & 1]->ch[c] = prv_node->ch[idx & 1]->ch[c];
            }
            now_node = now_node->ch[idx & 1];
            prv_node = prv_node->ch[idx & 1];
            idx >>= 1;
        }
    }

    T get(int time, int idx) {
        Node *node = roots[time];
        while (true) {
            if (idx == 0) return node->val;
            node = node->ch[idx & 1];
            idx >>= 1;
        }
    }
};

}  // namespace rklib

#endif  // RK_PERSISTENT_ARRAY_HPP