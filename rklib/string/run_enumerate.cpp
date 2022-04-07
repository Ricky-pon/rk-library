void enum_run_rec(int l, int r, string &s, vector<tuple<int, int, int>> &runs){
    if(r-l <= 1) return;
    int m = (l+r) / 2;
    enum_run_rec(l, m, s, runs);
    enum_run_rec(m, r, s, runs);
    
    auto f = [&](bool rev){
        string t = s.substr(l, r-l);
        if(rev){
            reverse(t.begin(), t.end());
            m = l+r-m;
        }

        int len = r-l, mid = m-l;
        string tl = t.substr(0, mid);
        reverse(tl.begin(), tl.end());
        string tr = t.substr(mid, len-mid) + t;
        auto zl = z_algorithm(tl), zr = z_algorithm(tr);
        zl.push_back(0);

        for(int k=1; mid-k>=0; ++k){
            int li = m-k-zl[k], ri = m+min(r-m, zr[len-k]);
            if(rev){
                swap(li, ri);
                li = l+r - li;
                ri = l+r - ri;
            }
            if(ri-li < 2*k) continue;
            if(li > 0 && s[li-1] == s[li-1+k]) continue;
            if(ri < (int)s.size() && s[ri] == s[ri-k]) continue;
            runs.push_back(make_tuple(li, ri, k));
        }
    };

    f(false);
    f(true);
}

vector<tuple<int, int, int>> enum_run(string &s){
    vector<tuple<int, int, int>> runs;
    enum_run_rec(0, s.size(), s, runs);
    sort(runs.begin(), runs.end());
    vector<tuple<int, int, int>> ret;
    for(int i=0; i<(int)runs.size(); ++i){
        if(i>0 && get<0>(runs[i]) == get<0>(runs[i-1]) && get<1>(runs[i]) == get<1>(runs[i-1])){
            continue;
        }
        auto [l, r, t] = runs[i];
        ret.push_back(make_tuple(t, l, r));
    }
    return ret;
}