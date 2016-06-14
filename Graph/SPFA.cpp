bool SPFA(int s) {
    memset (cnt, 0, sizeof (cnt));
    memset (vis, false, sizeof (vis));
    memset (d, INF, sizeof (d));
    cnt[s] = 1; d[s] = 0; vis[s] = true;
    std::queue<int> que;
    que.push (s);
    while (!que.empty ()) {
        int u = que.front ();
        que.pop ();
        vis[u] = false;
        for (int i=head[u]; ~i; i=edge[i].nex) {
            Edge &e = edge[i];
            if (d[e.v] > d[u] + e.w) {
                d[e.v] = d[u] + e.w;
                if (!vis[e.v]) {
                    vis[e.v] = true;
                    que.push (e.v);
                    if (++cnt[e.v] > n) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}
