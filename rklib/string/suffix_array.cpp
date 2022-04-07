void induced_sort(vector<int> &vec, int val_range, vector<int> &sa, vector<bool> &sl, vector<int> &lms_idx){
    vector<int> l(val_range, 0), r(val_range, 0);
    for(int c: vec){
        if(c+1 < val_range) ++l[c+1];
        ++r[c];
    }
    partial_sum(l.begin(), l.end(), l.begin());
    partial_sum(r.begin(), r.end(), r.begin());

    fill(sa.begin(), sa.end(), -1);

    for(int i=lms_idx.size()-1; i>=0; --i){
        sa[--r[vec[lms_idx[i]]]] = lms_idx[i];
    }

    for(int i: sa)if(i >= 1 && sl[i-1]){
        sa[l[vec[i-1]]++] = i-1;
    }

    fill(r.begin(), r.end(), 0);
    for(int c: vec) ++r[c];
    partial_sum(r.begin(), r.end(), r.begin());
    for(int k=sa.size()-1, i=sa[k]; k>=1; --k, i=sa[k])if(i >= 1 && !sl[i-1]){
        sa[--r[vec[i-1]]] = i-1;
    }
}

vector<int> sa_is(vector<int> &vec, int val_range){
    const int n = vec.size();
    vector<int> sa(n), lms_idx;
    vector<bool> sl(n);

    sl[n-1] = false;
    for(int i=n-2; i>=0; --i){
        sl[i] = (vec[i] > vec[i+1] || (vec[i] == vec[i+1] && sl[i+1]));
        if(sl[i] && !sl[i+1]) lms_idx.push_back(i+1);
    }
    reverse(lms_idx.begin(), lms_idx.end());

    induced_sort(vec, val_range, sa, sl, lms_idx);

    vector<int> new_lms_idx(lms_idx.size()), lms_vec(lms_idx.size());
    for(int i=0, k=0; i<n; ++i)if(!sl[sa[i]] && sa[i] >= 1 && sl[sa[i]-1]){
        new_lms_idx[k++] = sa[i];
    }

    int cur = 0;
    sa[n-1] = cur;
    for(size_t k=1; k<new_lms_idx.size(); ++k){
        int i = new_lms_idx[k-1], j = new_lms_idx[k];
        if(vec[i] != vec[j]){
            sa[j] = ++cur;
            continue;
        }
        bool flag = false;
        for(int a=i+1, b=j+1; ; ++a, ++b){
            if(vec[a] != vec[b]){
                flag = true;
                break;
            }
            if((!sl[a] && sl[a-1]) || (!sl[b] && sl[b-1])){
                flag = !((!sl[a] && sl[a-1]) && (!sl[b] && sl[b-1]));
                break;
            }
        }
        sa[j] = (flag ? ++cur : cur);
    }
    for(size_t i=0; i<lms_idx.size(); ++i){
        lms_vec[i] = sa[lms_idx[i]];
    }

    if(cur+1 < (int)lms_idx.size()){
        auto lms_sa = sa_is(lms_vec, cur+1);

        for(size_t i=0; i<lms_idx.size(); ++i){
            new_lms_idx[i] = lms_idx[lms_sa[i]];
        }
    }

    induced_sort(vec, val_range, sa, sl, new_lms_idx);

    return sa;
}

vector<int> suffix_array(string s){
    s += '$';
    vector<int> vec(s.size());
    for(int i=0; i<s.size(); ++i) vec[i] = s[i];
    auto sa = sa_is(vec, 128);
    sa.erase(sa.begin());
    return sa;
}