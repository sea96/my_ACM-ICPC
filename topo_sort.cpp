/*
    *return true 有环; false 无环
*/
bool topo_sort() {
    memset (deg, 0, sizeof (deg));
    for (int i=1; i<=n; ++i) {
        for (int j=0; j<edge[i].size (); ++j) {
            deg[edge[i][j]]++;
        }
    }
    int cnt = 0;
    std::queue<int> que;
    for (int i=1; i<=n; ++i) {
        if (!deg[i]) {
            que.push (i);
        }
    }
    while (!que.empty ()) {
        int u = que.front (); que.pop ();
        ans[++cnt] = u;
        for (int i=0; i<edge[u].size (); ++i) {
            int v = edge[u][i];
            if (!(--deg[v])) {
                que.push (v);
            }
        }
    }
    return cnt != n;
}
