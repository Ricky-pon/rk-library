vector<int> calc_right(vector<int> &sa, vector<int> &lcp) {
    int n = sa.size();

    vector<int> R(n);
    R[sa[n - 1]] = 0;

    stack<int> st;
    st.push(n - 1);

    for (int i = n - 2; i >= 0; --i) {
        R[sa[i]] = lcp[i];
        while (!st.empty() && sa[i] < sa[st.top()]) {
            R[sa[i]] = min(R[sa[i]], R[sa[st.top()]]);
            st.pop();
        }
        if (st.empty()) R[sa[i]] = 0;
        st.push(i);
    }

    return R;
}

vector<int> lz_parsing(string &s) {
    int n = s.size();
    auto sa = suffix_array(s);
    auto lcp = lcp_array(s, sa);

    auto R = calc_right(sa, lcp);
    reverse(sa.begin(), sa.end());
    reverse(lcp.begin(), lcp.end());
    auto L = calc_right(sa, lcp);

    vector<int> lzp;
    for (int i = 0; i < n;) {
        lzp.push_back(i);
        int j = i + max(L[i], R[i]);
        if (j == i) j = i + 1;
        i = j;
    }
    return lzp;
}