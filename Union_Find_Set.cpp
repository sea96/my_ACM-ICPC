/*
    并查集(Union_Find_Set)，核心函数：Union ()，Find ()；O(logn)
    Union (): 按秩合并两个不相同集合；
    Find (): 找到两个结点各自的root，顺便实现路径压缩(有递归和非递归版本)
    功能：合并集合，查询是否在同一集合；剩余未连通点数
*/
struct UF   {
    int rt[N], rk[N];
    void init(void) {
        memset (rt, -1, sizeof (rt));
        memset (rk, 0, sizeof (rk));
    }
    int Find(int x) {
        return (rt[x] == -1) ? x : rt[x] = Find (rt[x]);
    }
    void Union(int x, int y)    {
        x = Find (x), y = Find (y);
        if (x == y) return ;
        if (rk[x] > rk[y])   {
            rt[y] = x;  rk[x] += rk[y] + 1;
        }
        else    {
            rt[x] = y;  rk[y] += rk[x] + 1;
        }
    }
    bool same(int x, int y) {
        return Find (x) == Find (y);
    }
}uf;
 
//迭代形式的路径压缩
int Find(int x) {
    int p = x;
    while (rt[x] != -1) x = rt[x];
    while (p != x)  {
        int t = rt[p];  rt[p] = x;  p = t;
    }
    return x;
}
