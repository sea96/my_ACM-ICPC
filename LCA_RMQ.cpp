/*
    *LCA -> RMQ: LCA (u, v) = E[id[u] <= i <= id[v]中dep[i]最小的i];
    *RMQ预处理复杂度O(nlogn)，每次询问O (1)
    *(各数组大小:dp[N<<1][D], F[N<<1], dep[N<<1], id[N])
*/
void DFS(int u, int fa, int d, int &k) {
    id[u] = k; E[k] = u; dep[k++] = d;
    for (int i=head[u]; ~i; i=edge[i].nex) {
        int v = edge[i].v;
        if (v == fa) continue;
        DFS (v, u, d + 1, k);
        E[k] = u; dep[k++] = d;
    }
}
int min_dep(int i, int j) {
    return dep[i] < dep[j] ? i : j;
}
void init_RMQ(int k) {
    for (int i=0; i<k; ++i) dp[i][0] = i;
    for (int j=1; (1<<j)<=k; ++j) {
        for (int i=1; i+(1<<j)<k; ++i) {
            dp[i][j] = min_dep (dp[i][j-1], dp[i+(1<<(j-1))][j-1]);
        }
    }
}
int RMQ(int l, int r) {
    int k = 0; while (1<<(k+1) <= r - l + 1) k++;
    return min_dep (dp[l][k], dp[r-(1<<k)+1][k]);
}
int LCA(int u, int v) {
    u = id[u]; v = id[v];
    if (u > v) std::swap (u, v);
    return E[RMQ (u, v)];
}
