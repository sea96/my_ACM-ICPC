/*
    *Floyd_Warshall：任意两点的最短路，动态规划思想，复杂度O(V^3)
    *可求传递闭包
*/
void Floyd_Warshall() {
    for (int k=1; k<=n; ++k) {
        for (int i=1; i<=n; ++i) {
            for (int j=1; j<=n; ++j) {
                d[i][j] = std::min (d[i][j], d[i][k] + d[k][j]);
                //g[i][j] = (g[i][j] || (g[i][k] && g[k][j]));
            }
        }
    }
}
