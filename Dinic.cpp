/*
   *最大流之Dinic算法：不停地构造层次图，然后用阻塞流增广。
   *时间复杂度为 O(n^2*m)。
   *如果容量为1，复杂度为 O(min(n^(2/3), m^(1/2)*m)。
   *对于二分图最大匹配这样的特殊图，复杂度为 O(n^(1/2)*m)
*/
struct Max_Flow {
    struct Edge {
        int from, to, cap, flow;
    };
    std::vector<Edge> edges;
    std::vector<int> G[N]; //保存每个结点的弧在edges里的序号
    int level[N], cur[N]; //结点在层次图中的等级。结点当前弧下标。
    int n, m, s, t;
    //初始化顶点个数n，边的个数m在加边时统计，s和t分别为源点和汇点
    void init(int n) {
        this->n = n;
        for (int i=0; i<=n; ++i) {
            G[i].clear ();
        }
        edges.clear ();
    }
    void add_edge(int from, int to, int cap) {
        edges.push_back ((Edge) {from, to, cap, 0});
        edges.push_back ((Edge) {to, from, 0, 0});
        m = edges.size ();
        G[from].push_back (m - 2);
        G[to].push_back (m - 1);
    }
    bool BFS() {
        std::fill (level, level+1+n, -1);
        std::queue<int> que;
        level[s] = 0; que.push (s);
        while (!que.empty ()) {
            int u = que.front (); que.pop ();
            for (int i=0; i<G[u].size (); ++i) {
                Edge &e = edges[G[u][i]];
                if (level[e.to] == -1 && e.cap > e.flow) {
                    level[e.to] = level[u] + 1;
                    que.push (e.to);
                }
            }
        }
        return level[t] != -1;
    }
    int DFS(int u, int a) {
        if (u == t || a == 0) {
            return a; //a表示当前为止所有弧的最小残量
        }
        int flow = 0, f;
        for (int &i=cur[u]; i<G[u].size (); ++i) {
            Edge &e = edges[G[u][i]];
            if (level[u] + 1 == level[e.to]
            && (f = DFS (e.to, std::min (a, e.cap - e.flow))) > 0) {
                e.flow += f;
                edges[G[u][i]^1].flow -= f;
                flow += f; a -= f;
                if (a == 0) {
                    break;
                }
            }
        }
        return flow;
    }
    int Dinic(int s, int t) {
        this->s = s; this->t = t;
        int flow = 0;
        while (BFS ())  {
            std::fill (cur, cur+1+n, 0);
            flow += DFS (s, INF);
        }
        return flow;
    }
};
