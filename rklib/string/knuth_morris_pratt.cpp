template <class T>
struct KMP {
    T p;
    vector<int> fail;

    KMP(T p) : p(p) {
        fail.resize(p.size() + 1);
        fail[0] = -1;
        for (int i = -1, j = 0; j < (int)p.size();) {
            while (i >= 0 && p[i] != p[j]) i = fail[i];
            ++i;
            ++j;
            fail[j] = (j < (int)p.size() && p[i] == p[j] ? fail[i] : i);
        }
    }

    vector<int> match(T &t) {
        vector<int> pos;
        for (int i = 0, j = 0; j < (int)t.size();) {
            while (i >= 0 && p[i] != t[j]) i = fail[i];
            ++i;
            ++j;
            if (i == (int)p.size()) {
                pos.push_back(j - i);
                i = fail[i];
            }
        }
        return pos;
    }
};