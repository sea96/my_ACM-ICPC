//最大流
struct ISAP {
    static const int V = N;
    struct Edge {...};
    std::vector<Edge> edges;
    std::vector<int> G[V];  //保存每个结点的弧在edges里的序号
    int d[V], cur[V];  //结点在层次图中的等级。结点当前弧下标。
    int n, m, s, t;  //n：点个数，m：边条数（加边时统计），s：源点，t：汇点
    int pre[V], num[V];
    void init(int n) {...}
    void add_edge(int from, int to, int cap) {...}
    void BFS() {
        memset(d, INF, sizeof(d));
        std::queue<int> que;
        que.push(t);
        d[t] = 0;
        while (!que.empty()) {
            int u = que.front();
            que.pop();
            for (int i=0; i<G[u].size (); ++i) {
                Edge &e = edges[G[u][i]^1];
                if (d[e.from] == INF && e.cap > e.flow) {
                    d[e.from] = d[u] + 1;
                    que.push(e.from);
                }
            }
        }
    }
    int augment() {
        int u = t, aug = INF;
        while (u != s) {
            Edge &e = edges[pre[u]];
            aug = min(aug, e.cap - e.flow);
            u = e.from;
        }
        u = t;
        while (u != s) {
            edges[pre[u]].flow += aug;
            edges[pre[u]^1].flow -= aug;
            u = edges[pre[u]].from;
        }
        return aug;
    }
    int max_flow(int s, int t) {
        this->s = s; this->t = t;
        int flow = 0;
        BFS();
        memset(num, 0, sizeof(num));
        for (int i=0; i<n; ++i) if (d[i] != INF) num[d[i]]++;  //!!!
        int u = s;
        memset(cur, 0, sizeof(cur));
        while (d[s] < n) {
            if (u == t) {
                flow += augment();
                u = s;
            }
            bool advanced = false;
            for (int i=cur[u]; i<G[u].size(); ++i) {
                Edge &e = edges[G[u][i]];
                if (e.cap > e.flow && d[u] == d[e.to]+1) {
                    advanced = true;
                    pre[e.to] = G[u][i];
                    cur[u] = i;
                    u = e.to;
                    break;
                }
            }
            if (!advanced) {  //retreat
                int m = n - 1;
                for (int i=0; i<G[u].size(); ++i) {
                    Edge &e = edges[G[u][i]];
                    if (e.cap > e.flow) m = min(m, d[e.to]);
                }
                if (!(--num[d[u]])) break;  //gap优化
                num[d[u]=m+1]++;
                cur[u] = 0;
                if (u != s) u = edges[pre[u]].from;
            }
        }
        return flow;
    }
};
