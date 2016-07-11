int mod;  //全局变量，矩阵乘法取模
/*
    *矩阵快速幂处理线性递推关系f(n)=a1f(n-1)+a2f(n-2)+...+adf(n-d)
*/
struct Matrix {
    int row, col;
    int arr[N][N];
    Matrix(int r=0, int c=0) {
        row = r; col = c;
        memset (arr, 0, sizeof (arr));
    }
    Matrix operator * (const Matrix &B) {
        Matrix ret(row, B.col);
        for (int i=0; i<row; ++i) {
            for (int j=0; j<B.col; ++j) {
                for (int k=0; k<col; ++k) {
                    ret.arr[i][j] = (ret.arr[i][j] + (ll) arr[i][k] * B.arr[k][j]) % mod;
                }
            }
        }
        return ret;
    }
    void unit(int n) {
        row = col = n;
        for (int i=0; i<n; ++i) {
            arr[i][i] = 1;
        }
    }
};
Matrix operator ^ (Matrix X, int n) {
    Matrix ret; ret.unit (X.col);
    while (n) {
        if (n & 1) {
            ret = ret * X;
        }
        X = X * X;
        n >>= 1;
    }
    return ret;
}

int main() {
    //input a[i], f[i] (1<=i<=d)
    if (n <= d) {
        //f[n];
    } else {
        Matrix Fn(d+1, d+1), Fd(d+1, 1);
        for (int i=0; i<Fn.row-1; ++i) {
            Fn.arr[i][i+1] = 1;
        }
        for (int i=1; i<Fn.col; ++i) {
            Fn.arr[Fn.row-1][i] = a[d-i+1];
        }
        for (int i=0; i<Fd.row; ++i) {
            Fd.arr[i][0] = f[i];
        }
        Fn = Fn ^ (n - d);
        Fn = Fn * Fd;
        //f[n] = Fn.arr[d][0];
    }
    return 0;
}
