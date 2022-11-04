#ifndef RK_SUFFIX_AUTOMATON_HPP
#define RK_SUFFIX_AUTOMATON_HPP

#include <cassert>
#include <map>
#include <string>
#include <vector>

namespace rklib {

template <class T>
struct SuffixAutomaton {
   public:
    SuffixAutomaton() {}
    SuffixAutomaton(std::string &s) {
        for (auto c : s) add(c);
    }
    SuffixAutomaton(std::vector<T> &v) {
        for (auto x : v) add(x);
    }

    int max_len(int state) {
        assert(0 <= state && state < int(len.size()));
        return len[state];
    }

    int trans(int state, T c) {
        assert(0 <= state && state < int(nxt.size()));
        return nxt[state].find(c) == nxt[state].end() ? 0 : nxt[state][c];
    }

    int next_state(int state, T c) {
        assert(0 <= state && state < int(nxt.size()));
        while (state != 0 && trans(state, c) == 0) {
            state = link[state];
        }
        return trans(state, c);
    }

    int suffix_link(int state) {
        assert(0 <= state && state < int(link.size()));
        return link[state];
    }

    void add(T c) {
        int now = size;
        create_node();
        len[now] = len[last] + 1;
        int p;
        for (p = last; p >= 0 && nxt[p].find(c) == nxt[p].end(); p = link[p]) {
            nxt[p][c] = now;
        }
        if (p < 0) {
            link[now] = 0;
        } else {
            int q = nxt[p][c];
            if (len[q] == len[p] + 1) {
                link[now] = q;
            } else {
                int clone = size;
                create_node();
                len[clone] = len[p] + 1;
                link[clone] = link[q];
                nxt[clone] = nxt[q];
                for (; p >= 0 && nxt[p][c] == q; p = link[p]) {
                    nxt[p][c] = clone;
                }
                link[q] = link[now] = clone;
            }
        }
        last = now;
    }

    bool is_substring(std::vector<T> &p) {
        int cur = 0;
        for (auto c : p) {
            if (nxt[cur].find(c) == nxt[cur].end()) return false;
            cur = nxt[cur][c];
        }
        return true;
    }

    bool is_substring(std::string &p) {
        std::vector<int> v(p.size());
        for (size_t i = 0; i < p.size(); i++) {
            v[i] = p[i];
        }
        return is_substring(v);
    }

   private:
    std::vector<int> len = {0}, link = {-1};
    std::vector<std::map<T, int>> nxt = {{}};
    int size = 1, last = 0;

    void create_node() {
        ++size;
        len.resize(size);
        link.resize(size);
        nxt.resize(size);
    }
};

}  // namespace rklib

#endif  // RK_SUFFIX_AUTOMATON_HPP