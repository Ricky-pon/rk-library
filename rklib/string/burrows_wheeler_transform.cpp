struct BWT{
    int n, width;
    string bwt;
    vector<int> cnt_smaller, sa;
    vector<vector<int>> cnt;

    BWT(string s, int w): n(s.size()+1), width(w){
        sa = suffix_array(s);
        s += '$';
        sa.insert(sa.begin(), n-1);
        bwt.resize(n);
        for(int i=0; i<n; ++i){
            bwt[i] = (sa[i] == 0 ? s[n-1] : s[sa[i]-1]);
        }

        cnt_smaller.resize(128, 0);
        cnt.resize(128);
        for(int i=0; i<128; ++i) cnt[i].resize((n+1)/width+1, 0);
        vector<int> table(128, 0);
        for(int i=0; i<n; ++i){
            if(bwt[i]+1 < 128) ++cnt_smaller[bwt[i]+1];
            ++table[bwt[i]];
            if((i+1)%width == 0){
                for(int c=0; c<128; ++c) cnt[c][(i+1)/width] = table[c];
            }
        }
        partial_sum(cnt_smaller.begin(), cnt_smaller.end(), cnt_smaller.begin());
    }

    int get_cnt(int c, int k){
        int ret = cnt[c][k/width];
        for(int i=k/width*width; i<k; ++i) ret += (bwt[i] == c);
        return ret;
    }

    pair<int, int> fm_index(string &s){
        int l = 0, r = n;
        for(int i=s.size()-1; i>=0; --i){
            l = cnt_smaller[s[i]] + get_cnt(s[i], l);
            r = cnt_smaller[s[i]] + get_cnt(s[i], r);
        }
        return {l, r};
    }
};