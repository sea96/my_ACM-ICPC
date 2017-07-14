void GuassJordan(Matrix A, int n) {
    int i, j, k, r;
    for (i=0; i<n; ++i) {
        r = i;
        for (j=i+1; j<n; ++j) if (fabs(A[j][i]) > fabs(A[r][i])) r = j;
        //放弃这一行，直接处理下一行
        if (fabs (A[r][i]) < EPS) continue;
        if (r != i) for (j=0; j<=n; ++j) swap(A[r][j], A[i][j]);
        //与除了第i行外的其他行进行消元
        for (k=0; k<n; ++k) {
            if (k == i) continue;
            for (j=n; j>=i; --j) A[k][j] -= A[k][i] / A[i][i] * A[i][j];
        }
    }
}
