/*
 *矩阵快速幂处理线性递推关系f(n)=a1f(n-1)+a2f(n-2)+...+adf(n-d)
 */
typedef vector<int> Vec;  //矩阵的一行用vector表示
typedef vector<Vec> Mat;  //矩阵用多行vector表示

//Mat A(row, Vec(col)); 行长度：A.size(), 列长度：A[0].size()。
Mat matrix_mul(const Mat &A, const Mat &B) {
    Mat ret(A.size(), Vec(B[0].size()));
    for (int i=0; i<A.size(); ++i)
        for (int j=0; j<B[0].size(); ++j)
            for (int k=0; k<A[0].size(); ++k)
                ret[i][j] = (ret[i][j] + (ll)A[i][k]*B[k][j]%MOD) % MOD;  //可优化
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
    //input a[i], f[i] (1<=i<=d)
    if (n <= d) {
        //f[n];
    } else {
        Mat Fn(d+1, Vec(d+1)), Fd(d+1, Vec(1));
        for (int i=0; i<Fn.size()-1; ++i) {
            Fn[i][i+1] = 1;
        }
        for (int i=1; i<Fn[0].size(); ++i) {
            Fn[Fn.size()-1][i] = a[d-i+1];
        }
        for (int i=0; i<Fd.size(); ++i) {
            Fd[i][0] = f[i];
        }
        Fn = matrix_pow(Fn, n - d);
        Fn = matrix_mul(Fn, Fd);
        //f[n] = Fn[d][0];
    }
    return 0;
}
