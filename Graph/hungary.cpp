// 0-base O(VE)
struct Hungary {
    int n, m;
    vector<int> edge[N];
    int left[N], right[N];
    bool S[N], T[N];

    void init(int n, int m) {
        this->n = n;
        this->m = m;
        for (int i=0; i<n; ++i) edge[i].clear();
    }
    void addEdge(int u, int v) {
        edge[u].push_back(v);
    }
    bool match(int u) {
        S[u] = true;
        for (int v: edge[u]) if (!T[v]) {
            T[v] = true;
            if (left[v] == -1 || match(left[v])) {
                left[v] = u;
                right[u] = v;
                return true;
            }
        }
        return false;
    }
    int run() {  // 二分图匹配
        memset(left, -1, sizeof(left));
        memset(right, -1, sizeof(right));
        int ret = 0;
        for (int u=0; u<n; ++u) {
            memset(T, false, sizeof(T));
            memset(S, false, sizeof(S));
            if (match(u)) ret++;
        }
        return ret;
    }
    int minCover(vector<int> &X, vector<int> &Y) {  // 最小点覆盖
        int ret = run();
        memset(S, false, sizeof(S));
        memset(T, false, sizeof(T));
        for (int u=0; u<n; ++u) if (right[u] == -1) match(u);
        for (int u=0; u<n; ++u) if (!S[u]) X.push_back(u);  // S的未标记点
        for (int v=0; v<m; ++v) if (T[v]) Y.push_back(v);  // T的已标记点
        return ret;
    }
};
