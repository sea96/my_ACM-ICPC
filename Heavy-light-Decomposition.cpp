int sz[N], dep[N], rt[N][D], pos[N], belong[N];
void DFS1(int u, int fa) {
    sz[u] = 1; dep[u] = dep[fa] + 1;
    rt[u][0] = fa;
    for (int i=head[u]; ~i; i=edge[i].nex) {
        int v = edge[i].v;
        if (v == fa) {
            continue;
        }
        DFS1 (v, u);
        sz[u] += sz[v];
    }
}
void DFS2(int u, int fa, int chain) {
    int k = 0;
    pos[u] = ++loc;
    belong[u] = chain;
    for (int i=head[u]; ~i; i=edge[i].nex) {
        int v = edge[i].v;
        if (v == fa) {
            continue;
        }
        if (sz[v] > sz[k]) {
            k = v;
        }
    }
    if (k == 0) {
        return ;
    }
    DFS2 (k, u, chain);
    for (int i=head[u]; ~i; i=edge[i].nex) {
        int v = edge[i].v;
        if (v == fa || v == k) {
            continue;
        }
        DFS2 (v, u, v);
    }
}
