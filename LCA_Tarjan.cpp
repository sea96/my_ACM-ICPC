/*
    *LCA离线处理，Tarjan算法，复杂度O (N+Q)
    *对询问次序按深搜时遍历到的节点顺序进行重组，并查集找祖先
    *ans[i]表示第i个询问的LCA
*/
void DFS(int u) {
    rt[u] = u;
    for (int i=head[u]; ~i; i=edge[i].nex) {
        int v = edge[i].v;
        if (rt[v] == -1) {
            DFS (v);
            rt[v] = u;
        }
    }
    for (int i=headq[u]; ~i; i=query[i].nex) {
        int v = query[i].v;
        if (rt[v] != -1) {
            ans[query[i].id] = Find (v);
        }
    }
}
