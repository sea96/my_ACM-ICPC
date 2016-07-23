//单点更新，区间查询
struct ZKW {
    int sum[N<<2];
    int M;
    void build(int n) {
        for (M=1; M<n+2; M<<=1);
        for (int i=1; i<=n; ++i) {
            scanf ("%d", &sum[i+M]);
        }
        for (int i=M-1; i; --i) {
            sum[i] = sum[i<<1] + sum[i<<1|1];
        }
    }
    int query(int l, int r) {
        int ret = 0;
        for (l+=M-1, r+=M+1; l^r^1; l>>=1, r>>=1) {
            if (~l&1) ret += sum[l^1];
            if (r&1) ret += sum[r^1];
        }
        return ret;
    }
    void modify(int i, int v) {
        for (i+=M; i>=1; i>>=1) {
            sum[i] += v;
        }
    }
};
