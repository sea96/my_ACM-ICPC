/*
*BZOJ 1452 [JSOI2009]Count
*一个n*m的方格，初始时每个格子有一个权值，2种操作：
*1. 改变一个格子的权值；
*2. 求一个子矩阵某种特定权值出现的个数。
*/
struct FenwickTree2D {
    int c[105][N][N], n, m;
    void init(int n, int m) {
        memset (c, 0, sizeof (c));
        this->n = n;
        this->m = m;
    }
    // k=a[x][y], z=add
    void updata(int k, int x, int y, int z) {
        for (int i=x; i<=n; i+=i&(-i)) {
            for (int j=y; j<=m; j+=j&(-j)) {
                c[k][i][j] += z;
            }
        }
    }
    // query k(color)
    int query(int k, int x, int y) {
        int ret = 0;
        for (int i=x; i; i-=i&(-i)) {
            for (int j=y; j; j-=j&(-j)) {
                ret += c[k][i][j];
            }
        }
        return ret;
    }
    // count a
    int count(int x1, int x2, int y1, int y2, int a) {
        return query(a, x2, y2) - query(a, x2, y1-1) - query(a, x1-1, y2) + query(a, x1-1, y1-1);
    }
};
// 树状数组1D（间接区间修改，区间查询）
struct FenwickTree {
   // ...
}B, C;
int A[N];  // A[i] = A[i-1] + a[i];
void modify(int l, int r, int v) {  // [l, r] +v
	B.modify(l, v);
	B.modify(r + 1, -v);
	C.modify(l, v * l);
	C.modify(r + 1, -v*(r+1));
}
int query(int l, int r) {  // a[l]+a[l+1]+...+a[r]
   int ret = A[r] - A[l-1];
   ret += (r+1)*B.query(r) - l*B.query(l-1);
   ret -= C.query(r) - C.query(l-1);
   return ret;
}
