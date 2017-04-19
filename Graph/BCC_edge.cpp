void Tarjan(int u, int fa) {
    low[u] = dfn[u] = ++dfs_clock;
    sta[top++] = u;
    for (int i=head[u]; ~i; i=edges[i].nex) {
        int v = edges[i].v;
        if (v == fa) continue;
        if (!dfn[v]) Tarjan(v, u);
        low[u] = min(low[u], low[v]);  //合并写法
    }
    if (low[u] == dfn[u]) {  //u是双连通分量的第一个点
        for (; ;) {
            int x = sta[--top];
            low[x] = low[u];  //区分不同的双连通分量
            //low[u] != low[v] 的边就是桥
            if (x == u) break;
        }
    }
}

void find_bcc() {
    memset (dfn, 0, sizeof(dfn));
    dfs_clock = top = 0;
    for (int i=1; i<=n; ++i) Tarjan(1, 0);
}
