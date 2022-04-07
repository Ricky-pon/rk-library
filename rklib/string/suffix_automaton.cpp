struct SuffixAutomaton {
    vector<int> len = {0}, link = {-1};
    vector<map<char, int>> nxt = {{}};
    int sz = 1, last = 0;

    SuffixAutomaton() {}

    SuffixAutomaton(string &s) {
        for (auto c : s) add(c);
    }

    void create_node() {
        ++sz;
        len.resize(sz);
        link.resize(sz);
        nxt.resize(sz);
    }

    void add(char c) {
        int now = sz;
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
                int clone = sz;
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

    bool is_substring(string &p) {
        int cur = 0;
        for (auto c : p) {
            if (nxt[cur].find(c) == nxt[cur].end()) return false;
            cur = nxt[cur][c];
        }
        return true;
    }
};