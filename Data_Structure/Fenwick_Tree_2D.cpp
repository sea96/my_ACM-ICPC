/*
    *BZOJ 1452 [JSOI2009]Count
    *一个n*m的方格，初始时每个格子有一个权值，2种操作：
    *1. 改变一个格子的权值；
    *2. 求一个子矩阵某种特定权值出现的个数。
*/
//2D
struct Fenwick_Tree_2D {
    int c[105][N][N], n, m;
    void init(int n, int m) {
        memset (c, 0, sizeof (c));
        this->n = n;
        this->m = m;
    }
    //k=a[x][y], z=add
    void updata(int k, int x, int y, int z) {
        for (int i=x; i<=n; i+=i&(-i)) {
            for (int j=y; j<=m; j+=j&(-j)) {
                c[k][i][j] += z;
            }
        }
    }
    //query k(color)
    int query(int k, int x, int y) {
        int ret = 0;
        for (int i=x; i; i-=i&(-i)) {
            for (int j=y; j; j-=j&(-j)) {
                ret += c[k][i][j];
            }
        }
        return ret;
    }
    //count a
    int count(int x1, int x2, int y1, int y2, int a) {
        return query (a, x2, y2) - query (a, x2, y1-1) - query (a, x1-1, y2) + query (a, x1-1, y1-1);
    }
};
