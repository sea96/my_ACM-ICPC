struct ZKW {
    int M;
    int sum[N<<2];
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
}B, C;

int main() {
	B.modify (l, v);
            C.modify (l, v * l);
            if (r < n) {
                B.modify (r + 1, -v);
                C.modify (r + 1, -v*(r+1));
            }
	return 0;
}
    for (int i=1; i<=n; ++i) {
        scanf ("%d", a+i);
        A[i] = A[i-1] + a[i];
    }
    int l, r, v;
    char op[2];
    while (q--) {
        scanf ("%s%d%d", op, &l, &r);
        if (op[0] == 'Q') {
            ll ans = A[r] - A[l-1];
            ans += (r+1) * B.query (l, r) + (r-l+1) * B.query (1, l - 1);
            ans -= C.query (l, r);
            printf ("%I64d\n", ans);
        } else {
            scanf ("%d", &v);
            
        }
    }
