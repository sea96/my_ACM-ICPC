/*
    *点双连通分量，时间复杂度O(V+E)
*/
struct Edge {
    int u, v, nex;
}edges[M<<1];
int head[N], etot;
int dfn[N], low[N];
int sta[N], top;
vector<int> bcc[N];
int bcc_no[N];
bool iscut[N];
int dfs_clock, bcc_cnt;

void addEdge(int u, int v) {
    edges[etot] = Edge{u, v, head[u]};
    head[u] = etot++;
}
void initGraph() {
    memset(head, -1, sizeof(head));
    etot = 0;
}
int Tarjan(int u, int fa) {
    int lowu = dfn[u] = ++dfs_clock;
    int child = 0, v;
    for (int i=head[u]; ~i; i=edges[i].nex) if (!dfn[v=edges[i].v]) {
        sta[top++] = i;
        child++;
        int lowv = Tarjan(v, u);
        lowu = min(lowu, lowv);  // 用后代的low函数更新
        if (lowv >= dfn[u]) {
            iscut[u] = true;
            bcc[++bcc_cnt].clear();
            for (; ;) {
                Edge &e = edges[sta[--top]];
                if (bcc_no[e.u] != bcc_cnt) {
                    bcc[bcc_cnt].push_back(e.u);
                    bcc_no[e.u] = bcc_cnt;
                }
                if (bcc_no[e.v] != bcc_cnt) {
                    bcc[bcc_cnt].push_back(e.v);
                    bcc_no[e.v] = bcc_cnt;
                }
                if (e.u == u && e.v == v) break;
            }
        }
    } else if (dfn[v] < dfn[u] && v != fa) {
        sta[top++] = i;
        lowu = min(lowu, dfn[v]);  // 用反向边更新u的low函数
    }
    if (fa < 0 && child == 1) iscut[u] = false;
    return lowu;
}

void findBcc() {
    memset(dfn, 0, sizeof(dfn));
    memset(iscut, 0, sizeof(iscut));
    memset(bcc_no, 0, sizeof(bcc_no));
    dfs_clock = top = bcc_cnt = 0;
    for (int i=1; i<=n; ++i) {
        if (!dfn[i]) Tarjan(i, -1);
    }
}
