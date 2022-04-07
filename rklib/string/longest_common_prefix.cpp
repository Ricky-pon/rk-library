vector<int> lcp_array(string &s, vector<int> &sa) {
    int n = s.size();
    vector<int> rnk(n);
    for (int i = 0; i < n; ++i) rnk[sa[i]] = i;
    vector<int> lcp(n - 1);
    int h = 0;
    for (int i = 0; i < n; ++i) {
        if (h > 0) --h;
        if (rnk[i] == 0) continue;
        for (int j = sa[rnk[i] - 1]; j + h < n && i + h < n; ++h) {
            if (s[i + h] != s[j + h]) break;
        }
        lcp[rnk[i] - 1] = h;
    }
    return lcp;
}