/*
 *树链剖分（例题：BZOJ 1036 [ZJOI2008]树的统计Count）
 *I. 把结点u的权值改为t
 *II. 询问从点u到点v的路径上的节点的最大权值
 *III. 询问从点u到点v的路径上的节点的权值和 注意：从点u到点v的路径上的节点包括u和v本身
 */
int sz[N], son[N], fa[N], dep[N];
void DFS1(int u, int pa) {
    sz[u] = 1;
    dep[u] = dep[pa] + 1;
    fa[u] = pa;
    for (int i=0; i<edge[u].size(); ++i) {
        int v = edge[u][i];
        if (v == pa) continue;
        DFS1(v, u);
        if (sz[v] > sz[son[u]]) son[u] = v;
        sz[u] += sz[v];
    }
}

int top[N], dfn[N], dfs_clock;
void DFS2(int u, int t) {
    dfn[u] = ++dfs_clock;
    top[u] = t;
    if (son[u] != 0) DFS2(son[u], t);
    for (int i=0; i<edge[u].size(); ++i) {
        int v = edge[u][i];
        if (v == fa[u] || v == son[u]) continue;
        DFS2(v, v);
    }
}

int LCA(int u, int v) {
    while (top[u] != top[v]) {
        if (dep[top[u]] < dep[top[v]]) swap(u, v);
        u = fa[top[u]];
    }
    return dep[u] < dep[v] ? u : v;
}

//segment_tree/fenwick_tree部分

int get_sum(int u, int lca) {
    int ret = 0;
    while (top[u] != top[lca]) {
        ret += query_sum(dfn[top[u]], dfn[u], 1, n, 1);
        u = fa[top[u]];
    }
    ret += query_sum(dfn[lca], dfn[u], 1, n, 1);
    return ret;
}

int get_max(int u, int lca) {
    int ret = -INF;
    while (top[u] != top[lca]) {
        ret = max(ret, query_max(dfn[top[u]], dfn[u], 1, n, 1));
        u = fa[top[u]];
    }
    ret = max(ret, query_max(dfn[lca], dfn[u], 1, n, 1));
    return ret;
}

void prepare() {
    sz[0] = son[0] = dep[0] = 0;
    dfs_clock = 0;
    DFS1(1, 0);
    DFS2(1, 1);
    for (int i=1; i<=n; ++i) {
        updata(dfn[i], a[i], 1, n, 1);
    }
}
