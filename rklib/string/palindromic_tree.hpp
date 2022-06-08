#ifndef RK_PALINDROMIC_TREE_HPP
#define RK_PALINDROMIC_TREE_HPP

#include <map>
#include <string>
#include <vector>

namespace rklib {

struct PalindromicTree {
   public:
    const int rootm1 = 0, root0 = 1;
    std::vector<int> len = {-1, 0}, link = {rootm1, rootm1};
    std::vector<std::map<char, int>> nxt = {{}, {}};
    std::string s;

    PalindromicTree() {}
    PalindromicTree(std::string& t) {
        for (auto c : t) add(c);
    }

    void add(char c) {
        while (last != rootm1) {
            int idx = int(s.size()) - len[last] - 1;
            if (idx >= 0 && s[idx] == c) break;
            last = link[last];
        }
        if (nxt[last].find(c) != nxt[last].end()) {
            last = nxt[last][c];
            s += c;
            return;
        }

        create_node();
        int now = sz - 1;
        len[now] = len[last] + 2;

        int p = link[last];
        while (p != rootm1) {
            int idx = int(s.size()) - len[p] - 1;
            if (idx >= 0 && s[idx] == c) break;
            p = link[p];
        }
        link[now] = (nxt[p].find(c) == nxt[p].end() ? root0 : nxt[p][c]);

        nxt[last][c] = now;
        last = now;
        s += c;
    }

    std::vector<int> count(std::string& t) {
        std::vector<int> res(sz, 0);

        int now = rootm1;
        for (int i = 0; i < int(t.size()); ++i) {
            char c = t[i];
            if (nxt[rootm1].find(c) == nxt[rootm1].end()) now = rootm1;
            while (now != rootm1) {
                int idx = i - len[now] - 1;
                if (idx >= 0 && t[idx] == c &&
                    nxt[now].find(c) != nxt[now].end())
                    break;
                now = link[now];
            }
            now = (nxt[now].find(c) == nxt[now].end() ? root0 : nxt[now][c]);
            ++res[now];
        }

        for (int i = int(res.size()) - 1; i >= 1; i--) {
            res[link[i]] += res[i];
        }

        return res;
    }

   private:
    int sz = 2, last = root0;

    void create_node() {
        ++sz;
        len.resize(sz);
        link.resize(sz);
        nxt.resize(sz);
    }
};

}  // namespace rklib

#endif  // RK_PALINDROMIC_TREE_HPP