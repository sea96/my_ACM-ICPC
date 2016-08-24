/*
    *点双连通分量，时间复杂度O(V+E)
*/
int dfn[N];
std::vector<int> bcc[N];
bool is_cut[N];
int sta[N];
int dfs_clock, top, bcc_cnt;

std::vector<int> edge[N];
int n;

int Tarjan(int u, int fa) {
    int lowu = dfn[u] = ++dfs_clock;
    int child = 0;
    sta[top++] = u;
    for (int v: edge[u]) {
        if (!dfn[v]) {
            child++;
            int lowv = Tarjan (v, u);
            lowu = std::min (lowu, lowv);
            if (lowv >= dfn[u]) {
                is_cut[u] = true;
                bcc[++bcc_cnt].clear ();
                for (; ;) {
                    int x = sta[--top];
                    bcc[bcc_cnt].push_back (x);
                    if (x == v) break;
                }
                bcc[bcc_cnt].push_back (u);
            }
        } else if (dfn[v] < dfn[u] && v != fa) {
            lowu = std::min (lowu, dfn[v]);
        }
    }
    if (fa < 0 && child == 1) is_cut[u] = false;
    return lowu;
}

void find_bcc() {
    memset (dfn, 0, sizeof (dfn));
    memset (is_cut, false, sizeof (is_cut));
    dfs_clock = top = bcc_cnt = 0;
    for (int i=1; i<=n; ++i) {
        if (!dfn[i]) Tarjan (i, -1);
    }
}
