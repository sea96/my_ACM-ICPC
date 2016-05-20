/*
    *Kruskal：并查集实现，记录两点和距离，按距离升序排序，O (ElogE)
*/
struct Edge {
    int u, v, w;
    bool operator < (const Edge &rhs) const {
        return w < rhs.w;
    }
}edge[E];
sort (edge+1, edge+1+m);
if (!uf.same (x, y))    uf.Union (x, y), ans += w;
