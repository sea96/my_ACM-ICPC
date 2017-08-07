/*
    *LCA -> RMQ: LCA (u, v) = E[id[u] <= i <= id[v]中dep[i]最小的i];
    *RMQ预处理复杂度O(nlogn)，每次询问O (1)
    *(各数组大小:dp[N<<1][D], F[N<<1], dep[N<<1], id[N])
*/
void DFS(int u, int fa, int d, int &k) {
    id[u] = k; E[k] = u; dep[k++] = d;
    for (int i=head[u]; ~i; i=edge[i].nex) {
        int v = edge[i].v;
        if (v == fa) continue;
        DFS (v, u, d + 1, k);
        E[k] = u; dep[k++] = d;
    }
}
int min_dep(int i, int j) {
    return dep[i] < dep[j] ? i : j;
}
void init_RMQ(int k) {
    for (int i=0; i<k; ++i) dp[i][0] = i;
    for (int j=1; (1<<j)<=k; ++j) {
        for (int i=1; i+(1<<j)<k; ++i) {
            dp[i][j] = min_dep (dp[i][j-1], dp[i+(1<<(j-1))][j-1]);
        }
    }
}
int RMQ(int l, int r) {
    int k = 0; while (1<<(k+1) <= r - l + 1) k++;
    return min_dep (dp[l][k], dp[r-(1<<k)+1][k]);
}
int LCA(int u, int v) {
    u = id[u]; v = id[v];
    if (u > v) std::swap (u, v);
    return E[RMQ (u, v)];
}

/*
    *LCA(倍增法，二分搜索)：rt[u][i](i<D, 2^(D-1)>N) 表示u的第2^i的祖先
    *LCA预处理复杂度O (logn)，每次询问O (logn)
    *DFS中要记录点的深度以及它的父亲，(dep[u] = d; rt[u][0] = fa;)
*/
void DFS(int u, int fa) {
    rt[u][0] = fa;
    for (int v: edge[u]) {
        if (v == fa) continue;
        dep[v] = dep[u] + 1;
        DFS(v, u);
    }
}

void init_LCA() {
    for (int j=1; j<D; ++j) {
        for (int i=1; i<=n; ++i) {
            rt[i][j] = rt[i][j-1] == -1 ? -1 : rt[rt[i][j-1]][j-1];
        }
    }
}

int LCA(int u, int v) {
    if (dep[u] < dep[v]) swap(u, v);
    for (int i=0; i<D; ++i) {
        if ((dep[u] - dep[v]) >> i & 1) u = rt[u][i];
    }
    if (u == v) return u;
    for (int i=D-1; i>=0; --i) {
        if (rt[u][i] != rt[v][i]) {
            u = rt[u][i];
            v = rt[v][i];
        }
    }
    return rt[u][0];
}

/*
    *LCA离线处理，Tarjan算法，复杂度O(N+Q)
    *对询问次序按深搜时遍历到的节点顺序进行重组，并查集找祖先
    *ans[i]表示第i个询问的LCA
*/
void DFS(int u) {
    rt[u] = u;
    for (int i=head[u]; ~i; i=edge[i].nex) {
        int v = edge[i].v;
        if (rt[v] == -1) {
            DFS (v);
            rt[v] = u;
        }
    }
    for (int i=headq[u]; ~i; i=query[i].nex) {
        int v = query[i].v;
        if (rt[v] != -1) {
            ans[query[i].id] = Find (v);
        }
    }
}
