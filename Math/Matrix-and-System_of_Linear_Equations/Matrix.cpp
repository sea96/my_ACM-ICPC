        memset (arr, 0, sizeof (arr));
    }
    void set_size(int row, int col) {
        this->row = row;
        this->col = col;
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
    void init(int n) {
        set_size (n, n);
        for (int i=0; i<n; ++i) {
            arr[i][i] = 1;
        }
    }    
};
Matrix operator ^ (Matrix X, int n) {
    Matrix ret; ret.init (X.col);
    while (n) {
        if (n & 1) {
            ret = ret * X;
        }
        X = X * X;
        n >>= 1;
    }
    return ret;
}
