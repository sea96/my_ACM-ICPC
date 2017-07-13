struct Edge {
    int u, v, nex;
}edges[M<<1];
int head[N], etot;
int dfn[N], dfs_clock;
int bcc_no[N], bcc_cnt;
int sta[N], top;
bool vis[M<<1];

int Tarjan(int u) {
    int lowu = dfn[u] = ++dfs_clock;
    sta[top++] = u;
    for (int i=head[u]; ~i; i=edges[i].nex) {
		int v = edges[i].v;
        if (!dfn[v]) {
        	vis[i] = vis[i^1] = true;
            lowu = min(lowu, Tarjan(v));
        } else if (!vis[i]) {  // 支持有重边的图
	        vis[i] = vis[i^1] = true;
        	lowu = min(lowu, dfn[v]);
        }
    }
    if (lowu == dfn[u]) {  //u是边双连通分量的第一个点
        ++bcc_cnt;
        int x;
        do {
            x = sta[--top];
			bcc_no[x] = bcc_cnt;
        } while (x != u);
    }
    return lowu;
}

void findBCC() {
    memset(dfn, 0, sizeof(dfn));
    memset(bcc_no, 0, sizeof(bcc_no));
    dfs_clock = bcc_cnt = top = 0;
    for (int i=1; i<=n; ++i) if (!dfn[i]) Tarjan(i);
}
