// dp[N][N][D][D]
void initRMQ(int n, int m) {
    for (int i=1; i<=n; ++i) for (int j=1; j<=m; ++j) dp[i][j][0][0] = a[i][j];
    for (int i=0; (1<<i)<=n; ++i) {
        for (int j=0; (1<<j)<=m; ++j) {
            if (i == 0 && j == 0) continue;
            for (int r=1; r+(1<<i)-1<=n; ++r) {
                for (int c=1; c+(1<<j)-1<=m; ++c) {
                    if (i == 0) dp[r][c][i][j] = max(dp[r][c][i][j-1], dp[r][c+(1<<(j-1))][i][j-1]);  //1D
                    else dp[r][c][i][j] = max(dp[r][c][i-1][j], dp[r+(1<<(i-1))][c][i-1][j]);
                }
            }
        }
    }
}
// (r1, c1)->(r2, c2)
int RMQ2D(int r1, int c1, int r2, int c2) {
    int kr = 0, kc = 0;
    while ((1<<(kr+1)) <= r2-r1+1) kr++;
    while ((1<<(kc+1)) <= c2-c1+1) kc++;
    int max1 = dp[r1][c1][kr][kc];  //↖
    int max2 = dp[r1][c2-(1<<kc)+1][kr][kc];  //↗
    int max3 = dp[r2-(1<<kr)+1][c1][kr][kc];  //↙
    int max4 = dp[r2-(1<<kr)+1][c2-(1<<kc)+1][kr][kc];  //↘
    return max(max(max1, max2), max(max3, max4));
}
