#include <bits/stdc++.h>
using namespace std;
  
const int N = 1e5 + 5;
const int D = 21;
const int INF = 0x3f3f3f3f;
int a[N<<1];
int n, m, q;

inline void Min(int &a, int b) { if (a > b) a = b; }

vector<int> edge[N<<1];
int sz[N<<1], dep[N<<1], belong[N<<1], dfn[N<<1], num[N<<1], loc;
int anc[N<<1][D];
void DFS2(int u, int fa, int chain) {
    dfn[u] = ++loc; num[loc] = u; belong[u] = chain;
    int k = 0;
    for (int v: edge[u]) if (v != fa) if (sz[v] > sz[k]) k = v;
    if (!k) return ;
    DFS2(k, u, chain);
    for (int v: edge[u]) if (!(v == fa || v == k)) DFS2(v, u, v);
}
void DFS(int u, int fa) {
    sz[u] = 1; dep[u] = dep[fa] + 1;
    anc[u][0] = fa;
    for (int i=1; i<D; ++i) anc[u][i] = anc[anc[u][i-1]][i-1];
    for (int v: edge[u]) if (v != fa) {
        DFS(v, u);
        sz[u] += sz[v];
    }
}

int LCA(int u, int v) {
    if (dep[u] < dep[v]) swap(u, v);
    for (int i=0; i<D; ++i) if ((dep[u]-dep[v]) >> i & 1) u = anc[u][i];
    if (u == v) return u;
    for (int i=D-1; i>=0; --i) if (anc[u][i] != anc[v][i]) {
        u = anc[u][i]; v = anc[v][i];
    }
    return anc[u][0];
}

//边双连通分量，时间复杂度O(V+E)
struct Edge {
    int u, v, nex;
}edges[N<<1];
bool vis[N<<1];
int head[N], etot;
int bcc_no[N];
int sta[N], top;
int dfs_clock, bcc_cnt;
void addEdge(int u, int v) {
    edges[etot] = Edge{u, v, head[u]};
    head[u] = ++etot;
}
int Tarjan(int u, int fa) {
    int lowu = dfn[u] = ++dfs_clock;
    int child = 0, v;
    for (int i=head[u]; i; i=edges[i].nex) if (!dfn[v=edges[i].v]) {
        sta[top++] = i;
        child++;
        int lowv = Tarjan(v, u);
        lowu = min(lowu, lowv);  // 用后代的low函数更新
        if (lowv >= dfn[u]) {
            ++bcc_cnt;
            for (; ;) {
                Edge &e = edges[sta[--top]];
                if (bcc_no[e.u] != bcc_cnt) {
                    edge[bcc_cnt].push_back(e.u);
                    edge[e.u].push_back(bcc_cnt);
                    bcc_no[e.u] = bcc_cnt;
                }
                if (bcc_no[e.v] != bcc_cnt) {
                    edge[bcc_cnt].push_back(e.v);
                    edge[e.v].push_back(bcc_cnt);
                    bcc_no[e.v] = bcc_cnt;
                }
                if (e.u == u && e.v == v) break;
            }
        }
    } else if (dfn[v] < dfn[u] && v != fa) {
        sta[top++] = i;
        lowu = min(lowu, dfn[v]);  // 用反向边更新u的low函数
    }
//    if (fa < 0 && child == 1) iscut[u] = false;
    return lowu;
}

void findBCC() {
    bcc_cnt = n;
    for (int i=1; i<=n; ++i) if (!dfn[i]) Tarjan(i, -1);
}
 
#define lc o << 1
#define rc o << 1 | 1
#define mid (l + r >> 1)
int mn[(N*2)<<2];
void build(int o, int l, int r) {
    if (l == r) { mn[o] = a[num[l]]; return ; }
    build(lc, l, mid); build(rc, mid+1, r);
    mn[o] = min(mn[o<<1], mn[o<<1|1]);
}
int pos, val, ql, qr;
void modify(int o, int l, int r) {
    if (l == r) { mn[o] = val; return ; }
    if (pos <= mid) modify(lc, l, mid);
    else modify(rc, mid+1, r);
    mn[o] = min(mn[o<<1], mn[o<<1|1]);
}
int query(int o, int l, int r) {
    if (ql > r || qr < l) return INF;
    if (ql <= l && r <= qr) return mn[o];
    return min(query(lc, l, mid), query(rc, mid+1, r));
}

int nn;
void prepare() {
    findBCC();
    nn = n + bcc_cnt;
    for (int i=n+1; i<=nn; ++i) a[i] = INF;
    sz[0] = dep[0] = 0;
    DFS(1, 0);
    loc = 0;
    DFS2(1, 0, 1);
    build(1, 1, nn);
}

int getMin(int u, int lca) {
    int ret = INF;
    while (belong[u] != belong[lca]) {
        if (dep[belong[u]] < dep[belong[lca]]) swap(u, lca);
        ql = dfn[belong[u]]; qr = dfn[u];
        ret = min(ret, query(1, 1, nn));
        u = anc[belong[u]][0];
    }
    ql = dfn[lca]; qr = dfn[u];
    if (dep[belong[u]] < dep[belong[lca]]) swap(u, lca);
    Min(ret, query(1, 1, nn));
    if (u > n) Min(ret, a[anc[u][0]]);  // 特判
    return ret;
}

int main() {
    scanf("%d%d%d", &n, &m, &q);
    for (int i=1; i<=n; ++i) scanf("%d", &a[i]);
    int u, v;
    for (int i=1; i<=m; ++i) {
        scanf("%d%d", &u, &v);
        addEdge(u, v);
        addEdge(v, u);
    }
    prepare();
    char op[5];
    for (int i=1; i<=q; ++i) {
        scanf("%s%d%d", op, &u, &v);
        if (op[0] == 'C') {
            pos = dfn[u]; val = v;
            modify(1, 1, nn);
        } else {
            int lca = LCA(u, v);
            printf("%d\n", min(getMin(u, lca), getMin(v, lca)));
        }
    }
    return 0;
}
//http://blog.csdn.net/FromATP/article/details/69389705
