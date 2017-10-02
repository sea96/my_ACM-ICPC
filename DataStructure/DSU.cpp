// 递归形式的路径压缩
int Find(int x) {
    return rt[x] == -1 ? x : rt[x] = Find(rt[x]);
}
// 迭代形式的路径压缩
int Find(int x) {
    int p = x;
    while (rt[x] != -1) x = rt[x];
    while (p != x) { int t = rt[p]; rt[p] = x;  p = t; }
    return x;
}
