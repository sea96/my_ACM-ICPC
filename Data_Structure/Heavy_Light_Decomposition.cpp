/*
    *树链剖分（例题：BZOJ 1036 [ZJOI2008]树的统计Count）
    *I. 把结点u的权值改为t
　　*II. 询问从点u到点v的路径上的节点的最大权值
　　*III. 询问从点u到点v的路径上的节点的权值和 注意：从点u到点v的路径上的节点包括u和v本身
*/
int sz[N], dep[N], rt[N][D], pos[N], belong[N];
void DFS1(int u, int fa) {
    sz[u] = 1; dep[u] = dep[fa] + 1;
    rt[u][0] = fa;
    for (int i=head[u]; ~i; i=edge[i].nex) {
        int v = edge[i].v;
        if (v == fa)
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
//LCA部分
//segment_tree/fenwick_tree部分
int get_sum(int u, int lca) {
    int ret = 0;
    while (belong[u] != belong[lca]) {
        ret += query_sum (pos[belong[u]], pos[u], 1, n, 1);
        u = rt[belong[u]][0];
    }
    ret += query_sum (pos[lca], pos[u], 1, n, 1);
    return ret;
}
int get_max(int u, int lca) {
    int ret = -INF;
    while (belong[u] != belong[lca]) {
        ret = std::max (ret, query_max (pos[belong[u]], pos[u], 1, n, 1));
        u = rt[belong[u]][0];
    }
    ret = std::max (ret, query_max (pos[lca], pos[u], 1, n, 1));
    return ret;
}
