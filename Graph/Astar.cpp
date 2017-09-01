struct Edge {
    int v, w, nex;
}edges[M];
int head[N], etot;
void addEdge(int u, int v, int w) {
    edges[++etot] = Edge{v, w, head[u]};
    head[u] = etot;
}

int d[N];
bool done[N];
void Dijkstra(int s) {
    memset(done, false, sizeof(done[0])*(n+1));
    memset(d, INF, sizeof(d[0])*(n+1));
    priority_queue<pii> pque;
    pque.push(pii(-(d[s]=0), s));
    while (!pque.empty()) {
        int u = pque.top().second;
        pque.pop();
        if (done[u]) continue;
        done[u] = true;
        for (int i=head[u]; i; i=edges[i].nex) {
            Edge &e = edges[i];
            if (d[e.v] > d[u] + e.w) {
                d[e.v] = d[u] + e.w;
                pque.push(pii(-d[e.v], e.v));
            }
        }
    }
}

// A*求第k短路，先在反图以终点为起点跑一遍最短路
// 次短路只要在Dijkstra算法稍加修改（二维数组）
struct Node {
    int g, v;
    bool operator < (const Node &rhs) const {
        return !(g+d[v] < rhs.g+d[rhs.v]);  // f = g + d[v]
    }
};
int Astar(int s, int t, int k) {
    if (d[s] == INF) return -1;
    int cnt = 0;
    priority_queue<Node> pque;
    pque.push(Node{0, s});
    while (!pque.empty()) {
        Node it = pque.top();
        pque.pop();
        if (it.v == t && ++cnt == k) return it.g;
        for (int i=head[it.v]; i; i=edges[i].nex) {
            pque.push(Node{it.g+edges[i].w, edges[i].v});
        }
    }
    return -1;
}
