//UVA 11297
struct Seg_Tree_2D {
    int mx[N<<2][N<<2], mn[N<<2][N<<2];
    int n, m;
    int xo, x1, y1, x2, y2, x, y, v;
    int maxv, minv;
    int xleaf;

    void query1D(int l, int r, int o) {
        if (y1 <= l && r <= y2) {
            maxv = std::max (maxv, mx[xo][o]);
            minv = std::min (minv, mn[xo][o]);
        } else {
            int mid = l + r >> 1;
            if (y1 <= mid) query1D (lson);
            if (y2 > mid) query1D (rson);
        }
    }

    void query2D(int l, int r, int o) {
        if (x1 <= l && r <= x2) {
            xo = o;
            query1D (1, m, 1);
        } else {
            int mid = l + r >> 1;
            if (x1 <= mid) query2D (lson);
            if (x2 > mid) query2D (rson);
        }
    }

    void modify1D(int l, int r, int o) {
        if (l == r) {
            if (xleaf) {
                mx[xo][o] = mn[xo][o] = v;
                return ;
            }
            //modify (x1-x2, y)
            mx[xo][o] = std::max (mx[xo<<1][o], mx[xo<<1|1][o]);
            mn[xo][o] = std::min (mn[xo<<1][o], mn[xo<<1|1][o]);
        } else {
            int mid = l + r >> 1;
            if (y <= mid) modify1D (lson);
            else modify1D (rson);
            //modify (x, y1-y2)
            mx[xo][o] = std::max (mx[xo][o<<1], mx[xo][o<<1|1]);
            mn[xo][o] = std::min (mn[xo][o<<1], mn[xo][o<<1|1]);
        }
    }

    void modify2D(int l, int r, int o) {
        if (l == r) {
            xo = o; xleaf = 1;
            modify1D (1, m, 1);
        } else {
            int mid = l + r >> 1;
            if (x <= mid) modify2D (lson);
            else modify2D (rson);
            xo = o; xleaf = 0;
            modify1D (1, m, 1);
        }
    }

    void query() {
        maxv = -INF; minv = INF;
        query2D (1, n, 1);
    }

    void modify() {
        modify2D (1, n, 1);
    }
};

Seg_Tree2D t;

int main() {
    scanf ("%d", &n);
    t.n = t.m = n;
    for (int i=1; i<=n; ++i) {
        for (int j=1; j<=n; ++j) {
            scanf ("%d", &t.v);
            t.x = i; t.y = j;
            t.modify ();
        }
    }
    char op[5];
    int q;
    scanf ("%d", &q);
    while (q--) {
        scanf ("%s", &op);
        if (op[0] == 'q') {
            scanf ("%d%d%d%d", &t.x1, &t.y1, &t.x2, &t.y2);
            t.query ();
            printf ("%d %d\n", t.maxv, t.minv);
        } else {
            scanf ("%d%d%d", &t.x, &t.y, &t.v);
            t.modify ();
        }
    }
    return 0;
}
