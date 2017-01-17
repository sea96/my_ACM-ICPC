/*
 *矩阵快速幂处理线性递推关系：f(n)=a1f(n-1)+a2f(n-2)+...+adf(n-d)+c（C是常数）
 */
typedef vector<ll> Vec;  //矩阵的一行用vector表示
typedef vector<Vec> Mat;  //矩阵用多行vector表示

//Mat A(row, Vec(col)); 行长度：A.size(), 列长度：A[0].size()。
Mat matrix_mul(const Mat &A, const Mat &B) {
    Mat ret(A.size(), Vec(B[0].size()));
    for (int i=0; i<A.size(); ++i)
        for (int j=0; j<A[0].size(); ++j) if (A[i][j])
            for (int k=0; k<B[0].size(); ++k) if (B[j][k])
                add_mod(ret[i][k], A[i][j]*B[j][k]%mod);
    return ret;
}

Mat matrix_pow(Mat X, int n) {
    Mat ret(X.size(), Vec(X.size()));
    for (int i=0; i<X.size(); ++i) ret[i][i] = 1;
    for (; n; n>>=1) {
        if (n & 1) ret = matrix_mul(ret, X);
        X = matrix_mul(X, X);
    }
    return ret;
}

int main() {
    //input a[i], f[i] (1<=i<=d), 
    Mat Fd(1, Vec(d+1));
    Mat X(d+1, Vec(d+1));
    for (int i=0; i<d; ++i) Fd[0][i] = f[d-i];
    for (int i=0; i<d; ++i) X[i][0] = a[i+1];
    for (int i=0; i<d-1; ++i) X[i][i+1] = 1;
    X[d][0] = c; X[d][d] = 1;
    Mat Fn = matrix_mul(Fd, matrix_pow(X, n-d));
    //f[n] = Fn[0][0];
    return 0;
}
