/*
    *LCA(倍增法，二分搜索)：rt[u][i](i<D=20) 表示u的第2^i的祖先
    *LCA预处理复杂度O (logn)，每次询问O (logn)
    *DFS中要记录点的深度以及它的父亲，(dep[u] = d; rt[u][0] = fa;)
*/
void DFS(int u, int fa, int k) {
    rt[u][0] = fa; dep[u] = k;
    for (auto v: G[u]) {
        if (v == fa) continue;
        DFS (v, u, k + 1);
    }
}

void init_LCA() {
    for (int j=1; j<D; ++j) {
        for (int i=1; i<=n; ++i) {
            rt[i][j] = rt[i][j-1] == 0 ? 0 : rt[rt[i][j-1]][j-1];
        }
    }
}

int LCA(int u, int v) {
    if (dep[u] < dep[v]) {
        std::swap (u, v);
    }
    for (int i=0; i<D; ++i) {
        if ((dep[u] - dep[v]) >> i & 1) {
            u = rt[u][i];
        }
    }
    if (u == v) {
        return u;
    }
    for (int i=D-1; i>=0; --i) {
        if (rt[u][i] != rt[v][i]) {
            u = rt[u][i];
            v = rt[v][i];
        }
    }
    return rt[u][0];
}
