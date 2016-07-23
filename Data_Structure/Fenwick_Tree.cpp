//树状数组（间接区间修改，区间查询）
struct Fenwick_Tree {
    int C[N];
    int n;
    void init(int n) {
        this->n = n;
        memset (C, 0, sizeof (C));
    }
    void modify(int i, int v) {
        for (; i<=n; i+=i&-i) C[i] += v;
    }
    int query(int i) {
        int ret = 0;
        for (; i>0; i-=i&-i) ret += C[i];
        return ret;
    }
}B, C;

int main() {
	//modify [l, r] +v
	B.modify (l, v);
	B.modify (r + 1, -v);
	C.modify (l, v * l);
	C.modify (r + 1, -v*(r+1));
    
    //query [l, r]
    int ans = A[r] - A[l-1];  //A[i] = A[i-1] + a[i];
    ans += (r + 1) * B.query (r) - l * B.query (l - 1);
    ans -= C.query (r) - C.query (l - 1);
	return 0;
}
