struct Node {
    int lc, rc;
    int s;
};
Node nd[N*20];
int root[N];
int nsz;

int new_node() {
    return nsz++;
}

void build(int &o, int l, int r) {
    o = new_node();
    nd[o] = (Node){0, 0, 0};
    if (l == r) return ;
    int mid = l + r >> 1;
    build(nd[o].lc, l, mid);
    build(nd[o].rc, mid+1, r);
    nd[o].s = nd[nd[o].lc].s + nd[nd[o].rc].s;
}

void modify(int &o, int past, int l, int r, int p) {
    o = new_node();
    nd[o] = nd[past];
    if (l == r) {
        nd[o].s += 1;
        return ;
    }
    int mid = l + r >> 1;
    if (p <= mid) modify(nd[o].lc, nd[past].lc, l, mid, p);
    else modify(nd[o].rc, nd[past].rc, mid+1, r, p);
    nd[o].s = nd[nd[o].lc].s + nd[nd[o].rc].s;
}

int query(int o, int past, int l, int r, int k) {
    if (l == r) return l;
    int mid = l + r >> 1;
    int sum = nd[nd[o].lc].s - nd[nd[past].lc].s;
    if (k <= sum) return query(nd[o].lc, nd[past].lc, l, mid, k);
    else return query(nd[o].rc, nd[past].rc, mid+1, r, k-sum);
}

//query [l, r] kth small number
void solve() {
    std::sort(sa+1, sa+1+n);  //sa[]: sorted a[]
    int nn = std::unique (sa+1, sa+1+n) - sa - 1;

    nsz = 0;
    build(root[0], 1, nn);
    
    for (int i=1; i<=n; ++i) {
        int pos = std::lower_bound(sa+1, sa+1+nn, a[i]) - sa;
        modify(root[i], root[i-1], 1, nn, pos);
    }
    int l, r, k;
    while (m--) {
        scanf("%d%d%d", &l, &r, &k);
        int id = query(root[r], root[l-1], 1, nn, k);
        printf("%d\n", sa[id]);
    }
}
