//O(n^3)
struct KM {
    int n;
    std::vector<int> edge[N];
    int w[N][N];
    int lx[N], ly[N];
    int left[N];
    bool S[N], T[N];
    int slack[N];

    void init(int n) {
        this->n = n;
        for (int i=0; i<n; ++i) edge[i].clear();
        memset(w, 0, sizeof(w));
    }

    void add_edge(int u, int v, int cost) {
        edge[u].push_back(v);
        w[u][v] = cost;
    }

    bool match(int u) {
        S[u] = true;
        for (int v: edge[u]) if (!T[v]) {
            int a = lx[u] + ly[v] - w[u][v];
            if (!a) {
                T[v] = true;
                if (left[v] == -1 || match(left[v])) {
                    left[v] = u;
                    return true;
                }
            } else if (slack[v] > a) slack[v] = a;
        }
        return false;
    }

    void updata() {
        int a = INF;
        for (int v=0; v<n; ++v)
            if (!T[v] && a > slack[v]) a = slack[v];
        for (int i=0; i<n; ++i) {
            if (S[i]) lx[i] -= a;

            if (T[i]) ly[i] += a;
            else slack[i] -= a;
        }
    }

    void solve() {
        memset(left, -1, sizeof(left));
        for (int i=0; i<n; ++i) {
            lx[i] = *std::max_element(w[i], w[i]+n);
            ly[i] = 0;
        }
        for (int u=0; u<n; ++u) {
            for (int v=0; v<n; ++v) slack[v] = INF;
            for (; ;) {
                memset(S, false, sizeof(S));
                memset(T, false, sizeof(T));
                if (match(u)) break;
                else updata ();
            }
        }
    }
};
