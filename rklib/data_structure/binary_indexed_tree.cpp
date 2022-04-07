template<typename T> struct BinaryIndexedTree{
    vector<T> node;

    BinaryIndexedTree(int n){
        node.resize(n+1, {});
    }

    void update(int i, T x){
        ++i;
        while(i < (int)node.size()){
            node[i] += x;
            i += (i & -i);
        }
    }

    T query(int i){
        ++i;
        lint ret = 0;
        while(i){
            ret += node[i];
            i -= (i & -i);
        }
        return ret;
    }

    T query(int l, int r){
        return query(r-1) - query(l-1);
    }
};