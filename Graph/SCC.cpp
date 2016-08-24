/*
    *强连通分量，时间复杂度O(V+E) 
*/
int dfn[N], scc_no[N];
std::vector<int> scc[N];
int sta[N];
int dfs_clock, top, scc_cnt;

std::vector<int> edge[N];
int n, m;

int Tarjan(int u) {
    int lowu = dfn[u] = ++dfs_clock;
    sta[top++] = u;
    for (int v: edge[u]) {
        if (!dfn[v]) {
            int lowv = Tarjan (v);
            lowu = std::min (lowu, lowv);
        } else if (!scc_no[v]) {
            lowu = std::min (lowu, dfn[v]);
        }
    }
    if (lowu == dfn[u]) {
        scc[++scc_cnt].clear ();
        for (; ;) {
            int x = sta[--top];
            scc[scc_cnt].push_back (x);
            scc_no[x] = scc_cnt;
            if (x == u) break;
        }
    }
    return lowu;
}

void find_scc() {
    memset (dfn, 0, sizeof (dfn));
    memset (scc_no, 0, sizeof (scc_no));
    dfs_clock = top = scc_cnt = 0;
    for (int i=1; i<=n; ++i) {
        if (!dfn[i]) Tarjan (i);
    }
}
