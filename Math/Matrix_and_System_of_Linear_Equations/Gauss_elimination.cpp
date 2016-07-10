typedef double Matrix[N][N];

//要求系数矩阵可逆
//A是增广矩阵，A[i][n]是常数bi
//运行结束后A[i]][n]对应第i个未知数的值
void Guass_elimination(Matrix A, int n) {
    int i, j, k, r;
    //消元过程
    for (i=0; i<n; ++i) {
        //选择一行与第i行交换
        r = i;
        for (j=i+1; j<n; ++j) {
            if (fabs (A[j][i]) > fabs (A[r][i])) {
                r = j;
            }
        }
        if (r != i) {
            for (j=0; j<=n; ++j) {
                std::swap (A[r][j], A[i][j]);
            }
        }
        //与第i+1~n行进行逆序枚举消元
        for (j=n; j>=i; --j) {
            for (k=i+1; k<n; ++k) {
                A[k][j] -= A[k][i] / A[i][i] * A[i][j];
            }
        }
    }
    //回代过程
    for (i=n-1; i>=0; --i) {
        for (j=i+1; j<n; ++j) {
            A[i][n] -= A[j][n] * A[i][j];
        }
        A[i][n] /= A[i][i];
    }
}
