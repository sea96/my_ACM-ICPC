/*
 * 最小费用最大流算法：每次在残余网络中沿着最短路增广(流最大前提下费用最小)
*/
struct MinCostMaxFlow {
    struct Edge {
        int from, to, cap, flow, cost;
    };
    vector<Edge> edges;
    vector<int> G[N];
    bool vis[N];
    int d[N], p[N], a[N]; //最短费用，上一条弧的序号，可改进量
    int n, m, s, t;
    
    void init(int n) {
        this->n = n;
        for (int i=0; i<=n; ++i) G[i].clear();
        edges.clear();
    }
    void addEdge(int from, int to, int cap, int cost) {
        edges.push_back(Edge{from, to, cap, 0, cost});
        edges.push_back(Edge{to, from, 0, 0, -cost});
        m = edges.size();
        G[from].push_back(m-2);
        G[to].push_back(m-1);
    }
    bool SPFA(int &flow, int &cost) {
        memset(d, INF, sizeof(d));
        memset(vis, false, sizeof(vis));
        memset(p, -1, sizeof(p));
        d[s] = 0; vis[s] = true; p[s] = 0; a[s] = INF;
        queue<int> que; que.push(s);
        while (!que.empty()) {
            int u = que.front(); que.pop();
            vis[u] = false;
            for (int i=0; i<G[u].size(); ++i) {
                Edge &e = edges[G[u][i]];
                if (e.cap > e.flow && d[e.to] > d[u] + e.cost) {
                    d[e.to] = d[u] + e.cost;
                    p[e.to] = G[u][i];
                    a[e.to] = min(a[u], e.cap - e.flow);
                    if (!vis[e.to]) {
                        vis[e.to] = true;
                        que.push (e.to);
                    }
                }
            }
        }
        if (d[t] == INF) return false;
        flow += a[t];
        cost += d[t] * a[t];
        int u = t;
        while (u != s) {
            edges[p[u]].flow += a[t];
            edges[p[u]^1].flow -= a[t];
            u = edges[p[u]].from;
        }
        return true;
    }
    void run(int s, int t, int &flow, int &cost) {
        this->s = s; this->t = t;
        flow = cost = 0;
        while (SPFA(flow, cost));
    }
};
