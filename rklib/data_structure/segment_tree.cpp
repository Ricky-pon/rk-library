template<typename T, typename E, typename F, typename G>
struct SegTree{
    int sz = 1, seq_sz;
    T et;
    F f;
    G g;
    vector<T> node;

    SegTree(int sz_, T et, F f, G g): seq_sz(sz_), et(et), f(f), g(g){
        while(sz<sz_) sz <<= 1;
        node.resize(sz<<1, et);
    }

    void build(vector<T> &a){
        rep(i, a.size()) node[i+sz] = a[i];
        rFor(i, sz, 1) node[i] = f(node[i<<1], node[(i<<1)+1]);
    }

    void build(T x){
        rep(i, seq_sz) node[i+sz] = x;
        rFor(i, sz, 1) node[i] = f(node[i<<1], node[(i<<1)+1]);
    }

    void update(int i, E x){
        i += sz;
        node[i] = g(node[i], x);
        i >>= 1;
        while(i){
            node[i] = f(node[i<<1], node[(i<<1)+1]);
            i >>= 1;
        }
    }

    T query(int l, int r){
        T vl = et, vr = et;
        for(l+=sz, r+=sz; l<r; l>>=1, r>>=1){
            if(l&1) vl = f(vl, node[l++]);
            if(r&1) vr = f(node[--r], vr);
        }
        return f(vl, vr);
    }

    int search_left(int l, const T val, function<bool(T, T)> cmp, function<bool(T, T)> check){
        T sum = et;
        for(l+=sz; ; l>>=1){
            if(check(f(sum, node[l]), val)){
                while(l < sz){
                    if(check(f(sum, node[l<<1]), val)) l <<= 1;
                    else{
                        sum = f(sum, node[l<<1]);
                        l = (l<<1) + 1;
                    }
                }
                return l-sz;
            }
            if(__builtin_popcount(l+1) == 1) return seq_sz;
            if(l&1) sum = f(sum, node[l++]);
        }
    }

    int search_right(int r, const T val, function<bool(T, T)> cmp, function<bool(T, T)> check){
        T sum = et;
        for(r+=sz; ; r>>=1){
            if(check(f(sum, node[r]), val)){
                while(r < sz){
                    if(check(f(node[(r<<1)+1], sum), val)) r = (r<<1) + 1;
                    else{
                        sum = f(node[(r<<1)+1], sum);
                        r<<=1;
                    }
                }
                return r-sz;
            }
            if(__builtin_popcount(r) == 1) return -1;
            if(!(r&1)) sum = f(node[r--], sum);
        }
    }

    int lower_bound_left(int l, const T val, function<bool(T, T)> cmp=less<>()){
        auto check = [&cmp](auto a, auto b){return !cmp(a, b);};
        return search_left(l, val, cmp, check);
    }

    int upper_bound_left(int l, const T val, function<bool(T, T)> cmp=less<>()){
        auto check = [&cmp](auto a, auto b){return cmp(b, a);};
        return search_left(l, val, cmp, check);
    }

    int lower_bound_right(int l, const T val, function<bool(T, T)> cmp=greater<>()){
        auto check = [&cmp](auto a, auto b){return !cmp(b, a);};
        return search_right(l, val, cmp, check);
    }

    int upper_bound_right(int l, const T val, function<bool(T, T)> cmp=greater<>()){
        auto check = [&cmp](auto a, auto b){return cmp(a, b);};
        return search_right(l, val, cmp, check);
    }
};