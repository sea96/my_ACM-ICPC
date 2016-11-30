/*
 *高斯消元，求解同余方程组
 *MOD为全局变量
 *运行结束后，A[i][m]对应第i个未知数的值
 */
int Guass_elimination(Mat &A) {
    int n = A.size(), m = A[0].size()-1, i, j, row, col;
    for (row=0, col=0; row<n && col<m; ++row, ++col) {
        //选择一行与row行交换
        for (i=row; i<n; ++i) if (A[i][col]) break;
        if (i == n) { row--; continue; }
        if (i != row) {
            for (j=0; j<=m; ++j) swap(A[row][j], A[i][j]);
        }
        //利用最小公倍数消元
        //浮点数直接逆向枚举消元，见训练指南代码
        for (i=row+1; i<n; ++i) {
            if (A[i][col]) {
                int lcm = LCM(A[row][col], A[i][col]);
                int t1 = lcm / A[i][col], t2 = lcm / A[row][col];
                for (j=col; j<=m; ++j) {
                    A[i][j] = ((A[i][j]*t1 - A[row][j]*t2) % MOD + MOD) % MOD;
                }
            }
        }
    }

    //无解
    for (i=row; i<n; ++i) if (A[i][m] != 0) return -1;
    //无穷多解
    if (row < m) return 1;
    //唯一解
    for (i=m-1; i>=0; --i) {
        int &ans = A[i][m];
        for (j=i+1; j<m; ++j) {
            if (A[i][j]) ans = (ans - A[i][j]*A[j][m]%MOD + MOD) % MOD;
        }
        //A[i][i]x=ans(mod MOD)
        int x, y, d;
        ex_GCD(A[i][i], MOD, x, y, d);  //d = GCD(A[i][i], MOD);
        ans = (x * (ans/d) % MOD + MOD) % MOD;  //ans % d != 0时，无解
    }
    return 0;
}
