/*
 *2SAT，暴力判断，注意从0开始
 */
struct Two_SAT {
    int n;
    std::vector<int> edge[N<<1];
    bool mark[N<<1];
    int sta[N<<1], top;

    void init(int n) {
        this->n = n;
        for (int i=0; i<=n*2; ++i) edge[i].clear ();
        memset (mark, false, sizeof (mark));
    }
    //x+valx = true or y+valy = true
    void add_edge(int x, int valx, int y, int valy) {
        x = x * 2 + valx;
        y = y * 2 + valy;
        edge[x^1].push_back (y);
        edge[y^1].push_back (x);
    }
    bool DFS(int x) {
        if (mark[x^1]) return false;
        if (mark[x]) return true;
        mark[x] = true;
        sta[top++] = x;
        for (int y: edge[x]) {
            if (!DFS (y)) return false;
        }
        return true;
    }
    bool judge() {
        for (int i=0; i<n*2; i+=2) {
            if (!mark[i] && !mark[i+1]) {
                top = 0;
                if (!DFS (i)) {
                    while (top > 0) mark[sta[--top]] = false;
                    if (!DFS (i+1)) return false;
                }
            }
        }
        return true;
    }
};
//还有强连通判断的算法
bool judge() {
    find_scc ();  //n*2
    for (int i=0; i<n; ++i) {
        if (scc_no[2*i] == scc_no[2*i+1]) return false;
    }
    return true;
}
