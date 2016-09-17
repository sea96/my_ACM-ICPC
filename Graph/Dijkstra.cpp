struct Dijkstra {
    struct Edge {
        int u, v, w;
    };
    int n, m;
    vector<Edge> edges;
    vector<int> id[N];
    int d[N], p[N];
    bool done[N];

    void init(int n) {
        this->n = n;
        edges.clear();
        for (int i=1; i<=n; ++i) id[i].clear();
    }
    void add_edge(int u, int v, int w) {
        edges.push_back((Edge){u, v, w});
        m = edges.size();
        id[u].push_back(m-1);
    }
    void dijkstra(int s) {
        memset(d, INF, sizeof(d));
        memset(done, false, sizeof(done));
        d[s] = 0;
        priority_queue<pair<int, int>> pque;
        pque.push(make_pair(-d[s], s));
        while (!pque.empty()) {
            int u = pque.top().second;
            pque.pop();
            if (done[u]) continue;
            done[u] = true;
            for (int i: id[u]) {
                Edge &e = edges[i];
                if (d[e.v] > d[u] + e.w) {
                    d[e.v] = d[u] + e.w;
                    p[e.v] = i;
                    pque.push(make_pair(-d[e.v], e.v));
                }
            }
        }
    }
};
