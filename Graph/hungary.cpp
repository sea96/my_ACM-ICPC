//O(VE)
bool DFS(int u) {
    for (auto v: edges[u]) {
        if (!vis[v]) {
            vis[v] = true;
            if (lk[v] == -1 || DFS (lk[v])) {
                lk[v] = u;
                return true;
            }
        }
    }
    return false;
}
int hungary() {
    int ret = 0;
    memset (lk, -1, sizeof (lk));
    for (int i=1; i<=n; ++i) {
        memset (vis, false, sizeof (vis));
        if (DFS (i)) ret++;
    }
    return ret;
}
