/*
    *后缀数组，倍增算法实现，复杂度O(nlogn)
    *sa[i]: 第i小的后缀是在字符串位置，即后缀sa[i]
    *rank[i]: 后追i在sa数组下标，即第rank[i]小
    *height[i]: LCP (suffix (sa[i-1], sa[i]))
*/
int sa[N], rank[N], height[N];
int ws[N], wa[N], wb[N];

bool cmp(int *r, int a, int b, int l) {
    return (r[a] == r[b] && r[a+l] == r[b+l]);
}
//r数组为读入的字符串，m = max (r[i]) + 1，一般字符128足够了
//n为strlen (s) + 1，加上最后一个'\0'
void DA(char *r, int n, int m = 128) {
    int i, j, p, *x = wa, *y = wb;
    for (i=0; i<m; ++i) ws[i] = 0;
    for (i=0; i<n; ++i) ws[x[i]=r[i]]++;
    for (i=1; i<m; ++i) ws[i] += ws[i-1];
    for (i=n-1; i>=0; --i) sa[--ws[x[i]]] = i;
    for (j=1, p=1; p<n; j<<=1, m=p) {
        for (p=0, i=n-j; i<n; ++i) y[p++] = i;
        for (i=0; i<n; ++i) if (sa[i] >= j) y[p++] = sa[i] - j;
        for (i=0; i<m; ++i) ws[i] = 0;
        for (i=0; i<n; ++i) ws[x[y[i]]]++;
        for (i=1; i<m; ++i) ws[i] += ws[i-1];
        for (i=n-1; i>=0; --i) sa[--ws[x[y[i]]]] = y[i];
        std::swap (x, y);
        for (p=1, x[sa[0]]=0, i=1; i<n; ++i) {
            x[sa[i]] = cmp (y, sa[i-1], sa[i], j) ? p - 1 : p++;
        }
    }
}
void calc_height(char *r, int *sa, int n) {
    int i, j, k = 0;
    for (i=1; i<=n; ++i) rank[sa[i]] = i; //sa[0] = n(s[n]='\0')
    for (i=0; i<n; ++i) { //i: 后缀i
        if (k) k--;
        j = sa[rank[i]-1];
        while (r[i+k] == r[j+k]) k++;
        height[rank[i]] = k;  //其实并没有计算height[n]
    }
}
