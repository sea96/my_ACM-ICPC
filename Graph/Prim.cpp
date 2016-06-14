/*
    *Prim：Dijkstra思想，优先队列优化，适合稀疏图，O (ElogV)
    *不连通返回-1，或返回最小生成树长度(MST)
*/
int Prim(int s) {
    memset (vis, false, sizeof (vis));
    memset (d, INF, sizeof (d));
    std::priority_queue<Edge> que;
    for (int i=head[s]; ~i; i=edge[i].nex) {
        int v = edge[i].v, w = edge[i].w;
        if (d[v] > w) {
            d[v] = w;
            que.push (Edge (v, d[v]));
        }
    }
    vis[s] = true; d[s] = 0;
    int ret = 0;
    while (!que.empty ()) {
        int u = que.top ().v; que.pop ();
        if (vis[u]) {
        	continue;
        }
        vis[u] = true;
        ret += d[u];
        for (int i=head[u]; ~i; i=edge[i].nex) {
            Edge &e = edge[i];
            if (!vis[e.v] && d[e.v] > e.w) {
                d[e.v] = e.w;
                que.push (Edge (e.v, d[e.v]));
            }
        }
    }
    return ret;
}
