void Dijkstra(int s) {
    memset (vis, false, sizeof (vis));
    memset (d, INF, sizeof (d));
    d[s] = 0;
    std::priority_queue<Edge> pque;
    pque.push (Edge (s, d[s]));
    while (!pque.empty ()) {
        int u = pque.top ().v; pque.pop ();
        if (vis[u]) {
            continue;
        }
        vis[u] = true;
        for (int i=head[u]; ~i; i=edge[i].nex) {
            Edge &e = edge[i];
            if (!vis[e.v] && d[e.v] > d[u] + e.w) {
                d[e.v] = d[u] + e.w;
                pque.push (Edge (e.v, d[e.v]));
            }
        }
    }
}
