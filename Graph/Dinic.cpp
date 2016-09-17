/*
   *最大流之Dinic算法：不停地构造层次图，然后用阻塞流增广。
   *时间复杂度为 O(n^2*m)。
   *如果容量为1，复杂度为O(min(n^(2/3), m^(1/2)*m)。
   *对于二分图最大匹配这样的特殊图，复杂度为O(n^(1/2)*m)
*/
struct Dinic {
    struct Edge {
        int from, to, cap, flow;
    };
    vector<Edge> edges;
    vector<int> id[N];  //每个节点的弧在edges里的下标
    int d[N], cur[N];  //层次图中，到起点的距离，当前弧下标
    int n, m, s, t;  //节点数，边数，源点编号，汇点编号

    void init(int n) {
        this->n = n;  //n个点，0~n-1，新增的点下标从n开始
        for (int i=0; i<n; ++i) {
            id[i].clear();
        }
        edges.clear();
    }
    void add_edge(int from, int to, int cap) {
        edges.push_back((Edge){from, to, cap, 0});
        edges.push_back((Edge){to, from, 0, 0});
        m = edges.size();
        id[from].push_back(m-2);
        id[to].push_back(m-1);
    }
    bool BFS() {
        memset(d, -1, sizeof(d));
        d[s] = 0;
        queue<int> que;
        que.push(s);
        while (!que.empty()) {
            int u = que.front();
            que.pop();
            for (int i: id[u]) {
                Edge &e = edges[i];
                if (d[e.to] == -1 && e.cap > e.flow) {
                    d[e.to] = d[u] + 1;
                    que.push(e.to);
                }
            }
        }
        return d[t] != -1;
    }
    int DFS(int u, int a) {
        if (u == t || a == 0) return a;  //当前为止所有弧的最小残量
        int flow = 0, f;
        for (int &i=cur[u]; i<id[u].size(); ++i) {
            Edge &e = edges[id[u][i]];
            if (d[u] + 1 == d[e.to] 
                    && (f = DFS(e.to, min(a, e.cap-e.flow))) > 0) {
                e.flow += f;
                edges[id[u][i]^1].flow -= f;
                flow += f;
                a -= f;
                if (a == 0) break;
            }
        }
        return flow;
    }
    int max_flow(int s, int t) {
        this->s = s; this->t = t;
        int flow = 0;
        while (BFS()) {
            memset(cur, 0, sizeof(cur));
            flow += DFS(s, INF);
        }
        return flow;
    }
};

