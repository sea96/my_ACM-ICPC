/****************************************************************************/
									目录
排序
	1. 快速排序
	2. 归并排序

字符串处理
	1. KMP

数学
	1. 欧几里德算法(GCD)
		1.1 最大公约数
		1.2 扩展欧几里德算法
		1.3 逆元素
		1.4 快速GCD
		1.5 模线性方程
	2. 素数
	3. 快速幂
	4. 欧拉函数
	5. 约瑟夫环公式
	6. 第K个全排列
	
DP
	1. 最长上升子序列
		1.1 LIS	(O(n^2))
		1.2 LIS	(O(nlogn))
	2. 最长公共子序列
	3. 最长公共上升子序列

数据结构
	1. 线段树
		1.1 单点更新
		1.2 成段更新
	2. 树状数组
	3. 并查集
图论
	1. 最短路
		1.1 Dijkstra(邻接矩阵)
		1.2 Dijkstra(优先队列)
		1.3 Bellman_Ford(邻接表)
		1.4 SPFA(队列)
		1.5 Floyd_Warshall
	2. 拓扑排序(队列)
	3. 最小生成树
		3.1 Kruskal(并查集)
		3.2 Prim
STL
	1. lower_bound/upper_bound
	2. next_permutation
高精度
	1. 完全大数模板 By Kuangbin
	2. 大数阶乘求位数
计算几何
	1. 基本函数
	2. 凸包
黑科技
	1. 快速读入输出(读入输出外挂)！
	2. 解决爆栈，手动加栈！
	3. #include <bits/stdc++.h>
/****************************************************************************/

/****************************************************************************/
排序

1.快速排序
/*
	QuickSort：快速排序，随机标杆，O(nlogn)
*/
int Partition(int *a, int p, int r)
{
	int k = rand () % (r - p);
	swap (a[p+k], a[r]);
	int x = a[r];	int i = p;
	for (int j=p; j<r; ++j)
	{
		if (a[j] <= x)
		{
			swap (a[i++], a[j]);
		}
	}
	swap (a[i], a[r]);

	return i;
}

void QuickSort(int *a, int p, int r)
{
	if (p < r)
	{
		int q = Partition (a, p, r);
		QuickSort (a, p, q-1);
		QuickSort (a, q+1, r);		
	}
}

2.归并排序
/*
	MergeSort：归并排序，O(nlogn)
*/
void Merge(int *a, int p, int q, int r)
{
	int n1 = q - p + 1;
	int n2 = r - q;
	int i, j;
	for (i=1; i<=n1; ++i)	L[i] = a[p+i-1];
	for (j=1; j<=n2; ++j)	R[j] = a[q+j];
	L[n1+1] = R[n2+1] = INF;

	i = j = 1;
	for (int k=p; k<=r; ++k)		//求逆序数：cnt += n1 - i + 1;
		a[k] = (L[i] < R[j]) ? L[i++] : R[j++];
}

void MergeSort(int *a, int p, int r)
{
	if (p < r)
	{
		int q = p + (r - p) / 2;
		MergeSort (a, p, q);
		MergeSort (a, q + 1, r);
		Merge (a, p, q, r);
	}
}
/****************************************************************************/

/****************************************************************************/
字符串处理
1. KMP
/* 
 * next[]的含义：x[i-next[i]...i-1]=x[0...next[i]-1] By Kuangbin
 * next[i]为满足x[i-z...i-1]=x[0...z-1]的最大z值（就是x的自身匹配） 
 */ 
void kmp_pre(char x[],int m,int next[]) 
{ 
  int i,j; 
  j=next[0]=-1; 
  i=0; 
  while(i<m) 
  { 
    while(-1!=j && x[i]!=x[j])j=next[j]; 
    next[++i]=++j; 
  } 
} 
/* 
 * kmpNext[]的意思：next'[i]=next[next[...[next[i]]]] (直到next'[i]<0或者
	x[next'[i]]!=x[i]) 
 * 这样的预处理可以快一些 
 */ 
void preKMP(char x[],int m,int kmpNext[]) 
{ 
  int i,j; 
  j=kmpNext[0]=-1; 
  i=0; 
  while(i<m) 
  { 
    while(-1!=j && x[i]!=x[j])j=kmpNext[j]; 
    if(x[++i]==x[++j])kmpNext[i]=kmpNext[j]; 
    else kmpNext[i]=j; 
  } 
}
/* 
 * 返回x在y中出现的次数，可以重叠 
 */ 
int next[10010]; 
int KMP_Count(char x[],int m,char y[],int n) 
{//x是模式串，y是主串 
  int i,j; 
  int ans=0; 
  //preKMP(x,m,next); 
  kmp_pre(x,m,next); 
  i=j=0; 
  while(i<n) 
  { 
    while(-1!=j && y[i]!=x[j])j=next[j]; 
    i++;j++; 
    if(j>=m) 
    { 
      ans++; 
      j=next[j]; 
    } 
  } 
  return ans; 
} 
/****************************************************************************/

/****************************************************************************/
数学

1. 欧几里德算法

1.1 最大公约数
/*
	辗转相除法(欧几里德算法)
*/
ll GCD(ll a, ll b)	{return (b == 0) ? a : GCD (b, a % b);}

1.2 扩展欧几里德算法
/*
	扩展欧几里得算法， 求x, y 使得 gcd(a, b) = a * x + b * y;
	返回d = gcd(a,b)，和对应于等式 ax + by = d 中的x，y
	By Kuangbin
*/
ll ext_GCD(ll a, ll b, ll& x, ll& y)
{
    if (a == 0 && b == 0)	return -1;	//无最大公约数
    if (b == 0)	{x = 1;	y = 0;	return a;}
    ll d = ext_GCD (b, a % b, y, x);
    y -= a / b * x;
    return d;
}

1.3 逆元素
/*
    求逆元素，ax = 1(mod n)	By Kuangbin
*/
ll mid_reverse(ll a, ll n)
{
    ll x, y;
    ll d = ext_GCD (a, n, x, y);
    if (d == 1)	return (x % n + n) % n;
    else	return -1;
}

1.4 快速GCD
/*
	快速GCD，By 吉林模板
*/
int kgcd(int a, int b)
{
	if (a == 0) return b;
	if (b == 0) return a;
	if (!(a & 1) && !(b & 1)) return kgcd(a>>1, b>>1) << 1;
	else if (!(b & 1)) return kgcd(a, b>>1);
	else if (!(a & 1)) return kgcd(a>>1, b);
	else return kgcd(abs(a - b), min(a, b));
｝

1.5 模线性方程
/*
	模线性方程 a * x = b (% n)	By 吉林模板
*/
void modeq(int a, int b, int n) // ! n > 0
{
    int e, i, d, x, y;
    d = extgcd(a, n, x, y);
    if (b % d > 0) printf("No answer!\n");
    else
    {
        e = (x * (b / d)) % n;
        for (i = 0; i < d; i++) // !!! here x maybe < 0
            printf("%d-th ans: %d\n", i+1, (e+i*(n/d))%n);
    }
}

2. 素数
/*
	素性测试，输入正数，复杂度O (sqrt(n))
	By TiaoZhan
*/
bool is_prime(int n)
{
    for (int i=2; i*i<=n; ++i)
	{
		if (n % i == 0)	return false;
	}
	return n != 1;	//1例外
}

/*
	约数枚举，复杂度O (sqrt(n))
	By TiaoZhan
*/
vector<int> divisor(int n)
{
    vector<int> res;
    for (int i=1; i*i<=n; ++i)
	{
		if (n % i == 0)
		{
			res.push_back (i);
			if (n / i != i)	res.push_back (n / i);
		}
	}

	return res;
}

/*
	整数分解，复杂度O(sqrt(n))
	By TiaoZhan
*/
map<int, int> prime_factor(int n)
{
    map<int, int> res;
    for (int i=2; i*i<=n; ++i)
	{
        while (n % i == 0)	{++res[i];	n /= i;}
	}
	if (n != 1)	res[n] = 1;

	return res;
}

/*
    埃氏筛法：返回n以内素数的个数, 复杂度O (nloglogn)
	By TiaoZhan
*/
int prime[MAXN];
bool is_prime[MAXN];

int seive(int n)
{
    int p = 0;
    for (int i=0; i<=n; ++i)	is_prime[i] = false;
    is_prime[0] = is_prime[1] = false;
    for (int i=2; i<=n; ++i)
	{
        if (is_prime[i])
		{
            prime[++p] = i;			//有误
			for (int j=2*i; j<=n; ++j)	is_prime[j] = false;
		}
	}

	return p;
}

3. 快速幂
/*
	快速幂，复杂度O(logn)	By TiaoZhan
*/
ll pow_mod(ll x, ll n, ll mod)	//版本1
{
    ll res = 1;
    while (n > 0)
	{
		if (n & 1)	res = res * x % mod;
		x = x * x % mod;
		n >>= 1;
	}

	return res;
}

ll pow_mod(ll x, ll n, ll mod)	//版本2
{
    if (n == 0)	return 1;
    ll res = mod_pow (x * x % mod, n / 2, mod);
    if (n & 1)	res = res * x % mod;

    return res;
}

4. 欧拉函数
/*
	对正整数n，欧拉函数是少于或等于n的数中与n互质的数的数目。
	例如euler(8)=4，因为1,3,5,7均和8互质。
    Euler函数表达通式：euler(x)=x(1-1/p1)(1-1/p2)(1-1/p3)(1-1/p4)…(1-1/pn),
	其中p1,p2……pn为x的所有素因数，x是不为0的整数。euler(1)=1 
    欧拉公式的延伸：一个数的所有质因子之和是euler(n)*n/2。
*/
/*
	递推求数 欧拉函数 phi(i)	By 吉林模板
*/
for (i = 1; i <= maxn; i++) phi[i] = i;
for (i = 2; i <= maxn; i += 2) phi[i] /= 2;
for (i = 3; i <= maxn; i += 2) if (phi[i] == i)
    {
        for (j = i; j <= maxn; j += i)
            phi[j] = phi[j] / i * (i - 1);
    }
/*
	单独求数 欧拉函数 phi(x)	By 吉林模板
*/
unsigned euler(unsigned x)
{
    //  就是公式
    unsigned i, res=x;
    for (i = 2; i < (int)sqrt(x * 1.0) + 1; i++)
        if (x%i==0)
        {
            res = res / i * (i - 1);
            while (x % i == 0) x /= i; //  保证i 一定是素数
        }
    if (x > 1) res = res / x * (x - 1);
    return res;
}

5. 约瑟夫环公式
/*
	经典的约瑟夫(Joseph)水题
	(转载)约瑟夫环的递推公式：
    1. f[1]=0;　f[i]=(f[i-1]+m)%i; (i>1)	(0~n-1)
    2. f[1]=1;　f[i]=(f[i-1]+m)%i  (i>1);   if(f[i]==0) f[i]=i;		(1~n)
    3. P(1, m, k)=1 (i = 1);   P(i, m, k)=[P(i - 1, m, k ) + m - 1] % i + 1； 
		(i > 1, 此处先减1是为了让模i的值不为0)	(1~n)
*/

6. 第K个全排列
/*
	题意：求第K个全排列
    组合数学：首先，使用next_permutation 函数会超时，思路应该转变，
    摘抄网上的解法如下：
        假设第一位是a，不论a是什么数，axxxxxxxx一共有8!种选择。 
        297192 div 8! = 7，余14952，所以第一位是1-9中的第8个数，也就是8。 
        14952 div 7! = 2，余4872，所以第二位是3。 
        4872 div 6! = 6，余552，所以是第三位是1245679中的第7个，也就是 9。 
        552 div 5! = 4，余72，所以是124567中的第5个，也就是6。 
        72 div 4! = 2，余24，所以是4。 
        这时候就不用算了，因为24 = 4!，而剩下的数就是1257这4个，
		组成的排列的第24个必然是7521。
    以上解法只符合没有重复的序列，但是思路一致，
	把除法改为减法，每一次更新之后的全排列的数量
    即 Ann / Amm 的个数，可以用DFS实现
*/
long long fact(int x)
{
    long long res = 1;
    for (int i=1; i<=x; ++i)    res *= i;

    return res;
}

long long f(int step)
{
    long long res = fact (step);
    for (int i=0; i<26; ++i)    if (cnt[i])    res /= fact (cnt[i]);

    return res;    
}

void DFS(int step, long long k)
{
    if (step == len)
    {
        for (int i=0; i<len; ++i)    printf ("%c", 'A' + pos[i]);
        puts ("");    return ;
    }

    for (int i=0; i<26; ++i)
    {
        if (cnt[i] == 0)    continue;
        cnt[i]--;
        long long tmp = f (len - step - 1);
        if (tmp < k)    {k -= tmp;    cnt[i]++;}
        else
        {
            pos[step] = i;
            DFS (step+1, k);
            return ;
        }
    }
}
/****************************************************************************/

/****************************************************************************/
DP

1. 最长上升子序列

1.1 LIS	(O(n^2))
/*
    LIS(Longest Increasing Subsequence)：
        状态转移方程：dp[i] = max (dp[i], dp[j] + 1)   (1 <= j < i)
        附带有print输出路径函数
*/
void LIS(void)
{
	for (int i=1; i<=n; ++i)
    {
        dp[i] = 1;  rt[i] = -1;
        for (int j=1; j<i; ++j)
        {
            if (a[j] < a[i])
            {
                if (dp[i] < dp[j] + 1)
                {
                    dp[i] = dp[j] + 1;
                    rt[i] = j;
                }
            }
        }
        res = max (res, dp[i]);
    }
}

void print(int id, int n)
{
    if (rt[id] != -1)	print (rt[id], n);
    printf ("%d%c", a[id], (id == n) ? '\n' : ' ');
}

1.2 LIS	(O(nlogn))
/*
	LIS: 二分查找优化, O(nlogn)
		设当前最长递增子序列为len，考虑元素a[i]; 
        若d[len]<a[i]，则len++，并将d[len]=a[i]，那么d[]是有序的; 
        否则,在d[1~len]中二分查找，找到第一个比它小的元素d[j]，并用a[i]替换
*/
int LIS(void)
{
    int len = 1;
    int j;
    d[1] = a[1];
    
    for (int i=2; i<=n; ++i)
    {
        if (d[len] < a[i])
            d[++len] = a[i];
        else
        {
            j = lower_bound (d+1, d+1+len, a[i]) - d;
            d[j] = a[i];
        }
    }
    return len;
}

2. 最长公共子序列
/*
    LCS(Longest Common Subsequence)：
        状态转移方程：dp[i][j] = dp[i-1][j-1] + 1; (s[i-1] == t[i-1])
					dp[i][j] = max (dp[i][j-1], dp[i-1][j]);(s[i-1] != t[i-1])
		附带有print输出路径函数
*/
void LCS(void)
{
    for (int i=1; i<=len_s; ++i)		//数组范围：0~len-1
    {
        for (int j=1; j<=len_t; ++j)
        {
            if (s[i-1] == t[j-1])
            {
                dp[i][j] = dp[i-1][j-1] + 1;	rt[i][j] = 0;
            }
            else if (dp[i-1][j] >= dp[i][j-1])
            {
                dp[i][j] = dp[i-1][j];	rt[i][j] = 1;
            }
            else
            {
                dp[i][j] = dp[i][j-1];	rt[i][j] = -1;
            }
        }
    }
}

void print(int x, int y)
{
    if (!x && !y)    return ;

    if (rt[x][y] == 0)
    {
        print (x-1, y-1);    printf ("%c", s[x-1]);
    }
    else if (rt[x][y] == 1)
    {
        print (x-1, y);    printf ("%c", s[x-1]);
    }
    else
    {
        print (x, y-1);    printf ("%c", t[y-1]);
    }
}

3. 最长公共上升子序列
/*
    LCIS：if(a[i]==b[j]) dp[j]=max(dp[j],dp[k]);  (1<=k<j && b[k]<b[j])
            打印路径时按照b[i]来输出
*/
void print(int x, int y, int last)
{
    if (x == 0 || y == 0)    return ;

    print (rx[x][y], ry[x][y], y);
    if (y != last)
    {
        ++cnt;
        if (cnt == 1)    printf ("%d", b[y]);
        else    printf (" %d", b[y]);
    }
}

void LCIS(int n, int m)
{
    memset (dp, 0, sizeof (dp));
    memset (rx, 0, sizeof (rx));
    memset (ry, 0, sizeof (ry));

    int ans = 0;    int k = 0;
    int sx = 0, sy = 0;
    for (int i=1; i<=n; ++i)
    {
        k = 0;
        for (int j=1; j<=m; ++j)
        {
            dp[i][j] = dp[i-1][j];
            rx[i][j] = i - 1;    ry[i][j] = j;
            if (a[i] == b[j])
            {
                if (dp[i][j] < dp[i][k] + 1)
                {
                    dp[i][j] = dp[i][k] + 1;
                    rx[i][j] = i;    ry[i][j] = k;
                }
            }
            else if (a[i] > b[j])
            {
                if (dp[i][j] > dp[i][k])    k = j;
            }
            if (ans < dp[i][j])
            {
                ans = dp[i][j];
                sx = i;    sy = j;
            }
        }
    }

    cnt = 0;
    printf ("%d\n", ans);
    print (sx, sy, -1);    puts ("");
}
/****************************************************************************/

/****************************************************************************/
数据结构

1. 线段树

1.1 单点更新
/*
	线段树：单点更新，结构体包含树的各种属性，复杂度O(logn)
	功能：区间查询，单点更新
	注意：结构体开2倍大小
*/
#define lson l, mid, rt << 1
#define rson mid + 1, r, rt << 1 | 1

const int MAXN = 1e5 + 10;
const int INF = 0x3f3f3f3f;
struct Node
{
	int v, mx, mn, sum, len;
}node[MAXN << 2];

void push_up(int rt)
{
	node[rt].sum = node[rt<<1].sum + node[rt<<1|1].sum;
}

void build(int l, int r, int rt)
{
	if (l == r)	{scanf ("%d", &node[rt].sum);	return ;}
	int mid = (l + r) >> 1;
	build (lson);
	build (rson);

	push_up (rt);
}

void updata(int p, int add, int l, int r, int rt)
{
	if (l == r)	{node[rt].sum += add;	return ;}
	int mid = (l + r) >> 1;
	if (p <= mid)	updata (p, add, lson);
	else	updata (p, add, rson);

	push_up (rt);
}

int query(int ql, int qr, int l, int r, int rt)
{
	if (ql <= l && r <= qr)	{return node[rt].sum;}
	int ans = 0;
	int mid = (l + r) >> 1;
	if (ql <= mid)	ans += query (ql, qr, lson);
	if (qr > mid)	ans += query (ql, qr, rson);

	return ans;
}

1.2 成段更新
/*
	线段树：成段更新，延迟标记，即每次更新不更新到底，延迟到下一次更新or查询
*/
#define lson l, mid, rt << 1
#define rson mid + 1, r, rt << 1 | 1

const int MAXN = 1e5 + 10;
const int INF = 0x3f3f3f3f;
struct Node
{
	int v, sum, add, mx, mn, len;
}node[MAXN<<2];

void push_up(int rt)
{
	node[rt].sum = node[rt<<1].sum + node[rt<<1|1].sum;
}

void push_down(int rt, int c)
{
	if (node[rt].add)
	{
		node[rt<<1].add += node[rt].add;
		node[rt<<1|1].add += node[rt].add;
		node[rt<<1].sum += node[rt].add * (c - (c >> 1));
		node[rt<<1|1].sum += node[rt].add * (c >> 1);
		node[rt].add = 0;
	}
}

void build(int l, int r, int rt)
{
	node[rt].add = 0;
	if (l == r)	{scanf ("%d", &node[rt].sum);	return ;}
	int mid = (l + r) >> 1;
	build (lson);	build (rson);

	push_up (rt);
}

void updata(int ql, int qr, int c, int l, int r, int rt)
{
	if (ql <= l && r <= qr){node[rt].sum += c*(r-l+1); node[rt].add += c; return ;}

	push_down (rt, r - l + 1);

	int mid = (l + r) >> 1;
	if (ql <= mid)	updata (ql, qr, c, lson);
	if (qr > mid)	updata (ql, qr, c, rson);

	push_up (rt);
}

int query(int ql, int qr, int l, int r, int rt)
{
	if (ql <= l && r <= qr)	{return node[rt].sum;}

	push_down (rt, r - l + 1);

	int mid = (l + r) >> 1;		int ans = 0;
	if (ql <= mid)	ans += query (ql, qr, lson);
	if (qr > mid)	ans += query (ql, qr, rson);

	return ans;
}

2. 树状数组
/*
	实现求[l, r]的和
	注意：初始化为0，c[i] 用add (i, v)添加并修改，
	[l, r]的和:sum (r) - sum (l - 1)
*/
int c[MAXN];

int low_bit(int x)	{return x & (-x);}

void add(int i, int v)
{
    while (i <= n)	{c[i] += v;	i += low_bit (i);}
}

int sum(int i)
{
    int res = 0;
    while (i > 0)	{res += c[i];	i -= low_bit (i);}
    return res;
}

3. 并查集
/*
	核心函数：Union ()，Find ()；复杂度O(logn)
	Union ()，合并两个不相同集合；
	Find ()，找到两个结点各自的root，顺便实现路径压缩
	功能：合并集合，查询是否在同一集合；剩余未连通点数
*/
void init(void)	{memset (rt, -1, sizeof (rt));}

int Find(int x)	{return (rt[x] == -1) ? x : rt[x] = Find (rt[x]);}

void Union(int x, int y)		//简化版
{
	x = Find (x);	y = Find (y);
	if (x > y)	rt[y] = x;
	else if (x < y)	rt[x] = y;
}

bool same(int x, int y)	{return (Find (x) == Find (y));}
/****************************************************************************/

/****************************************************************************/
图论

1. 最短路

1.1 Dijkstra(邻接矩阵)
/*
	Dijkstra：单源最短路，邻接矩阵实现，权值非负
	求出起点s到任意点的最短距离d[n]，O(n^2)
	顺带给出打印路径函数，rt[s] = -1;
*/
const int MAXN = 1e4 + 10;
const int INF = 0x3f3f3f3f;
int w[MAXN][MAXN];
int d[MAXN];
bool vis[MAXN];
int rt[MAXN];
int n, m;

void Dijkstra(int s)
{
	for (int i=1; i<=n; ++i)	d[i] = INF;		//初始化INF
	memset (vis, 0, sizeof (vis));
	memset (rt, -1, sizeof (rt));
	d[s] = 0;						//s到s的距离为0

	for (int i=1; i<=n; ++i)
	{
		int mn = INF;	int x = -1;
		for (int j=1; j<=n; ++j)
		{
			if (!vis[j] && d[j] < mn)	mn = d[x=j];
		}
		if (x == -1)	break;		//找出当前最合适的点，贪心思想
		vis[x] = 1;
		for (int j=1; j<=n; ++j)
		{
			if (d[j] > d[x] + w[x][j])		//松弛操作，更新root
			{
				d[j] = d[x] + w[x][j];	rt[j] = x;
			}
		}
	}
}

void print(int x)
{
	if (rt[x] != -1)	{print (rt[x]);	printf (" %d", x);}		//递归输出最短路上所有的点
	else	printf ("%d", x);		//这是第一个点
}

1.2 Dijkstra(优先队列)
/*
	Dijkstra：单源最短路，优先队列优化，O(ElogE)
	注意：对G[i].clear ();
*/
const int MAXN = 1e4 + 10;
const int INF = 0x3f3f3f3f;
struct Edge
{
	int v, w;
	Edge (int _v = 0, int _w = 0) : v (_v), w (_w) {}
	bool operator < (const Edge & r) const
	{
		return w > r.w;
	}
};
vector<Edge> G[MAXN];
bool vis[MAXN];
int d[MAXN];
int n, m;

void Dijkstra(int s)
{
	priority_queue<Edge> Q;
	memset (vis, 0, sizeof (vis));
	for (int i=1; i<=n; ++i)	d[i] = INF;
	d[s] = 0;	Q.push (Edge (s, 0));

	Edge x;
	while (!Q.empty ())
	{
		x = Q.top ();	Q.pop ();
		int u = x.v;
		if (vis[u])	continue;
		vis[u] = true;
		for (int i=0; i<G[u].size (); ++i)
		{
			int v = G[u][i].v;
			int w = G[u][i].w;
			if (!vis[v] && d[v] > d[u] + w)
			{
				d[v] = d[u] + w;
				Q.push (Edge (v, d[v]));
			}
		}
	}
}

void add_edge(int u, int v, int w)
{
	G[u].push_back (Edge (v, w));
}

1.3 Bellman_Ford(邻接表)
/*
	Bellman_Ford：单源最短路，邻接表实现，可处理负边权图
	return：true：有环；false：无环
	注意：对G.clear ();
*/
const int MAXN = 1e2 + 10;
const int MAXM = 2e4 + 10;
const int INF = 0x3f3f3f3f;
struct Edge
{
	int u, v, w;
	Edge (int _u = 0, int _v = 0, int _w = 0) : u (_u), v (_v), w (_w) {}
};
int d[MAXN];
int cnt[MAXN];
bool vis[MAXN];
vector<Edge> G;
int n, m;

bool Bellman_Ford(int s)
{
	for (int i=1; i<=n; ++i)	d[i] = INF;
	d[s] = 0;

	bool flag;
	for (int i=1; i<=n-1; ++i)		//最多更新n-1次
	{
		flag = false;
		for (int j=0; j<G.size (); ++j)
		{
			int u = G[j].u;	int v = G[j].v;	int w = G[j].w;
			if (d[v] > d[u] + w)
			{
				d[v] = d[u] + w;	flag = true;
			}
		}
		if (!flag)	break;		//若没更新，则直接结束
	}

	for (int i=1; i<=n-1; ++i)
	{
		for (int j=0; j<G.size (); ++j)		//若存在负环回路，则还会更新
		{
			int u = G[j].u;	int v = G[j].v;	int w = G[j].w;
			if (d[v] > d[u] + w)	return true;
		}
	}

	return false;
}

void add_edge(int u, int v, int w)
{
	G.push_back (Edge (u, v, w));
}

1.4 SPFA(队列)
/*
	SPFA:单源最短路，队列实现，复杂度不定，bellman_Ford 优化版
	return：true：有环；false：无环
*/
const int MAXN = 1e2 + 10;
const int MAXM = 2e4 + 10;
const int INF = 0x3f3f3f3f;
struct Edge
{
	int v, w;
	Edge (int _v = 0, int _w = 0) : v (_v), w (_w) {}
};
int d[MAXN];
int cnt[MAXN];
bool vis[MAXN];
vector<Edge> G[MAXN];
int n, m;

bool SPFA(int s)
{
	memset (cnt, 0, sizeof (cnt));
	memset (vis, false, sizeof (vis));
	for (int i=1; i<=n; ++i)	d[i] = INF;
	d[s] = 0;	cnt[s] = 1;	vis[s] = true;

	queue<int> Q;	Q.push (s);
	while (!Q.empty ())
	{
		int u = Q.front ();	Q.pop ();
		vis[u] = false;
		for (int i=0; i<G[u].size (); ++i)
		{
			int v = G[u][i].v;	int w = G[u][i].w;
			if (d[v] > d[u] + w)	d[v] = d[u] + w;		//有错
			if (!vis[v])
			{
				vis[v] = true;	Q.push (v);
				if (++cnt[v] > n)	return true;
			}
		}
	}

	return false;
}

void add_edge(int u, int v, int w)
{
	G[u].push_back (Edge (v, w));
}

1.5 Floyd_Warshall
/*
	Floyd_Warshall：任意两点的最短路，动态规划思想，复杂度O(n^3)
	可求传递闭包
*/
void Floyd_Warshall(void)
{
	for (int k=1; k<=n; ++i)
	{
		for (int i=1; i<=n; ++i)
		{
			for (int j=1; j<=n; ++j)
			{
				d[i][j] = min (d[i][j], d[i][k] + d[k][j]);
				//d[i][j] = (d[i][j] || (d[i][k] && d[k][j]));	//传递闭包
			}
		}
	}
}

2. 拓扑排序(队列)
/*
    TopoSort：拓扑排序
    in[]：每个点的入度；    ans[]：排序后的结果；
    return：true：有环；false：无环
*/
bool TopoSort(void)
{
    memset (in, 0, sizeof (in));        		//入度清空
    for (int i=1; i<=n; ++i)					//所有有指向的点的入度累加
        for (int j=0; j<G[i].size (); ++j)    in[G[i][j]]++;

    queue<int> Q;    int cnt = 0;
    for (int i=1; i<=n; ++i)    {if (!in[i]) Q.push (i);}	//入读为0的入队

    while (!Q.empty ())
    {
        int u = Q.front ();    Q.pop ();
        ans[++cnt] = u;							//ans[]存储拓扑排序后的点
        for (int j=0; j<G[u].size (); ++j)
        {
            int v = G[u][j];
            in[v]--;						//每次减少一个入度
            if (!in[v])    Q.push (v);		//若为0，应最先入队，拓扑序靠前
        }
    }

    if (cnt == n)    return false;
    else return true;                	//若没有n个点，表示有环，YES
}

3. 最小生成树

3.1 Kruskal(并查集)
/*
	Kruskal：并查集实现，建立结构体，记录两点和距离，依照距离升序排序
*/
bool cmp(Node x, Node y) {return x.w < y.w;}
sort (node+1, node+1+m, cmp);
if (!same (x, y))	{Union (x, y);	ans += w;}

3.2 Prim
/*
	Prim：Dijkstra实现，传入起点s，不连通返回-1，或返回最小生成树长度(MST)
*/
const int MAXN = 1e3 + 10;
const int INF = 0x3f3f3f3f;
int d[MAXN];
int w[MAXN][MAXN];
bool vis[MAXN];
int n, m;

int Prim(int s)
{
	memset (vis, false, sizeof (vis));		//有误
	for (int i=1; i<=n; ++i)	d[i] = INF;
	vis[s] = true;	d[s] = 0;	int ans = 0;

	for (int i=1; i<=n; ++i)
	{
		int mn = INF;	int x = -1;
        for (int j=1; j<=n; ++j)
		{
            if (!vis[j] && mn > d[j])	mn = d[x=j];
		}
		if (x == -1)	return -1;
		vis[x] = true;	ans += mn;
		for (int j=1;j<=n; ++j)
		{
            if (!vis[j] && d[j] > d[x] + w[x][j])
				d[j] = d[x] + w[x][j];
		}
	}

	return ans;
}
/****************************************************************************/



/****************************************************************************/
STL (Standard Template Library)
1. lower_bound/upper_bound
/*
	lower_bound (begin, end, key)：第一个大于或等于key的位置(指针)
	upper_bound (begin, end, key)：最后一个大于或等于key的位置
	注意：范围左闭右开，找不到返回end(即最后一个元素的后一个位置)
*/
int p1 = lower_bound (a+1, a+1+n, k) - a;
int p2 = upper_bound (a+1, a+1+n, k) - a;
int len = p2 - p1;		//相同元素的个数

2. next_permutation
/*
	数组改变成下一个排列的顺序
	注意：生成全排列要先 sort
*/
do{...}while (next_permutation (a+1, a+1+n));

3. *max_element (a+1, a+1+n)/ *min_element (a+1, a+1+n)
/****************************************************************************/

/****************************************************************************/
高精度
/* 
 * 完全大数模板 By Kuangbin
 * 输出cin>>a 
 * 输出a.print(); 
 * 注意：这个输入不能自动去掉前导0的，
 *		先读入到char数组，去掉前导0，再用构造函数。
 */ 
#define MAXN 9999 
#define MAXSIZE 1010 
#define DLEN 4 

class BigNum 
{ 
private: 
    int a[500];  //可以控制大数的位数 
    int len; 
public: 
    BigNum(){len=1;memset(a,0,sizeof(a));}  //构造函数 
    BigNum(const int);     //将一个int类型的变量转化成大数 
    BigNum(const char*);   //将一个字符串类型的变量转化为大数 
    BigNum(const BigNum &); //拷贝构造函数 
    BigNum &operator=(const BigNum &); //重载赋值运算符，大数之间进行赋值运算 
    friend istream& operator>>(istream&,BigNum&); //重载输入运算符 
    friend ostream& operator<<(ostream&,BigNum&); //重载输出运算符 
 
    BigNum operator+(const BigNum &)const;  //重载加法运算符，两个大数之间的相加运算 
    BigNum operator-(const BigNum &)const;  //重载减法运算符，两个大数之间的相减运算 
    BigNum operator*(const BigNum &)const;  //重载乘法运算符，两个大数之间的相乘运算 
    BigNum operator/(const int &)const;     //重载除法运算符，大数对一个整数进行相除
运算 
 
    BigNum operator^(const int &)const;     //大数的n次方运算 
    int operator%(const int &)const;        //大数对一个int类型的变量进行取模运算
    bool operator>(const BigNum &T)const;   //大数和另一个大数的大小比较 
    bool operator>(const int &t)const;      //大数和一个int类型的变量的大小比较
 
    void print();        //输出大数 
}; 
BigNum::BigNum(const int b)   //将一个int类型的变量转化为大数 
{ 
    int c,d=b; 
    len=0; 
    memset(a,0,sizeof(a)); 
    while(d>MAXN) 
    { 
        c=d-(d/(MAXN+1))*(MAXN+1); 
        d=d/(MAXN+1); 
        a[len++]=c; 
    } 
    a[len++]=d; 
} 
BigNum::BigNum(const char *s)  //将一个字符串类型的变量转化为大数 
{ 
    int t,k,index,L,i; 
    memset(a,0,sizeof(a)); 
    L=strlen(s); 
    len=L/DLEN; 
    if(L%DLEN)len++; 
    index=0; 
    for(i=L-1;i>=0;i-=DLEN) 
    { 
        t=0; 
        k=i-DLEN+1; 
        if(k<0)k=0; 
        for(int j=k;j<=i;j++) 
            t=t*10+s[j]-'0'; 
        a[index++]=t; 
    } 
} 
BigNum::BigNum(const BigNum &T):len(T.len)  //拷贝构造函数 
{ 
    int i; 
    memset(a,0,sizeof(a)); 
    for(i=0;i<len;i++) 
        a[i]=T.a[i]; 
} 
BigNum & BigNum::operator=(const BigNum &n)  //重载赋值运算符，大数之间赋值运算
{ 
    int i; 
    len=n.len; 
    memset(a,0,sizeof(a)); 
    for(i=0;i<len;i++) 
        a[i]=n.a[i]; 
    return *this; 
} 
istream& operator>>(istream &in,BigNum &b) 
{ 
    char ch[MAXSIZE*4]; 
    int i=-1; 
    in>>ch; 
    int L=strlen(ch); 
    int count=0,sum=0; 
    for(i=L-1;i>=0;) 
    { 
        sum=0; 
        int t=1; 
        for(int j=0;j<4&&i>=0;j++,i--,t*=10) 
        { 
            sum+=(ch[i]-'0')*t; 
        } 
        b.a[count]=sum; 
        count++; 
    } 
    b.len=count++; 
    return in; 
} 
ostream& operator<<(ostream& out,BigNum& b)  //重载输出运算符 
{ 
    int i; 
    cout<<b.a[b.len-1]; 
    for(i=b.len-2;i>=0;i--) 
    { 
        printf("%04d",b.a[i]); 
    } 
    return out; 
} 
BigNum BigNum::operator+(const BigNum &T)const   //两个大数之间的相加运
{ 
    BigNum t(*this); 
    int i,big; 
    big=T.len>len?T.len:len; 
    for(i=0;i<big;i++) 
    { 
        t.a[i]+=T.a[i]; 
        if(t.a[i]>MAXN) 
        { 
            t.a[i+1]++; 
            t.a[i]-=MAXN+1; 
        } 
    } 
    if(t.a[big]!=0) 
       t.len=big+1; 
    else t.len=big; 
    return t; 
} 
BigNum BigNum::operator-(const BigNum &T)const  //两个大数之间的相减运算 
{ 
    int i,j,big; 
    bool flag; 
    BigNum t1,t2; 
    if(*this>T) 
    { 
        t1=*this; 
        t2=T; 
        flag=0; 
    } 
    else 
    { 
        t1=T; 
        t2=*this; 
        flag=1; 
    } 
    big=t1.len; 
    for(i=0;i<big;i++) 
    { 
        if(t1.a[i]<t2.a[i]) 
        { 
            j=i+1; 
            while(t1.a[j]==0) 
                j++; 
            t1.a[j--]--; 
            while(j>i) 
                t1.a[j--]+=MAXN; 
            t1.a[i]+=MAXN+1-t2.a[i]; 
        } 
        else t1.a[i]-=t2.a[i]; 
    } 
    t1.len=big; 
    while(t1.a[len-1]==0 && t1.len>1) 
    { 
        t1.len--; 
        big--; 
    } 
    if(flag) 
        t1.a[big-1]=0-t1.a[big-1]; 
    return t1; 
} 
BigNum BigNum::operator*(const BigNum &T)const  //两个大数之间的相乘 
{ 
    BigNum ret; 
    int i,j,up; 
    int temp,temp1; 
    for(i=0;i<len;i++) 
    { 
        up=0; 
        for(j=0;j<T.len;j++) 
        { 
            temp=a[i]*T.a[j]+ret.a[i+j]+up; 
            if(temp>MAXN) 
            { 
                temp1=temp-temp/(MAXN+1)*(MAXN+1); 
                up=temp/(MAXN+1); 
                ret.a[i+j]=temp1; 
            } 
            else 
            { 
                up=0; 
                ret.a[i+j]=temp; 
            } 
        } 
        if(up!=0) 
           ret.a[i+j]=up; 
    } 
    ret.len=i+j; 
    while(ret.a[ret.len-1]==0 && ret.len>1)ret.len--; 
    return ret; 
} 
BigNum BigNum::operator/(const int &b)const  //大数对一个整数进行相除运算
{ 
    BigNum ret; 
    int i,down=0; 
    for(i=len-1;i>=0;i--) 
    { 
        ret.a[i]=(a[i]+down*(MAXN+1))/b; 
        down=a[i]+down*(MAXN+1)-ret.a[i]*b; 
    } 
    ret.len=len; 
    while(ret.a[ret.len-1]==0 && ret.len>1) 
        ret.len--; 
    return ret; 
} 
int BigNum::operator%(const int &b)const   //大数对一个 int类型的变量进行取模 
{ 
    int i,d=0; 
    for(i=len-1;i>=0;i--) 
        d=((d*(MAXN+1))%b+a[i])%b; 
    return d; 
} 
BigNum BigNum::operator^(const int &n)const  //大数的n次方运算 
{ 
    BigNum t,ret(1); 
    int i; 
    if(n<0)exit(-1); 
    if(n==0)return 1; 
    if(n==1)return *this; 
    int m=n; 
    while(m>1) 
    { 
        t=*this; 
        for(i=1;(i<<1)<=m;i<<=1) 
           t=t*t; 
        m-=i; 
        ret=ret*t; 
        if(m==1)ret=ret*(*this); 
    } 
    return ret; 
} 
bool BigNum::operator>(const BigNum &T)const    //大数和另一个大数的大小比较 
{ 
    int ln; 
    if(len>T.len)return true; 
    else if(len==T.len) 
    { 
        ln=len-1; 
        while(a[ln]==T.a[ln]&&ln>=0) 
          ln--; 
        if(ln>=0 && a[ln]>T.a[ln]) 
           return true; 
        else 
           return false; 
    } 
    else 
       return false; 
} 

bool BigNum::operator>(const int &t)const  //大数和一个int类型的变量的大小比较
{ 
    BigNum b(t); 
    return *this>b; 
} 
void BigNum::print()   //输出大数 
{ 
    int i; 
    printf("%d",a[len-1]); 
    for(i=len-2;i>=0;i--) 
      printf("%04d",a[i]); 
    printf("\n"); 
} 
BigNum f[110];//卡特兰数 
 
int main() 
{ 
  f[0]=1; 
  for(int i=1;i<=100;i++) 
    f[i]=f[i-1]*(4*i-2)/(i+1);//卡特兰数递推式 
  int n; 
  while(scanf("%d",&n)==1) 
  { 
    if(n==-1)break; 
    f[n].print(); 
  } 
  return 0; 
} 

2. 大数介乘求位数
/*
	log10(n!) = log10(1*2*..*n) = log10(1) + log10(2) + ...+log10(n)
	解释：123456=1.23456*10^5;
		  log10(123456)=5.09151;
		  log10(1.23456*10^5)=log10(1.23456)+log10(10^5)=0.09151+5;
		  故int(log10(n))+1 就是n的位数
*/
void work(int n)
{
	double ans = 0;
	for (int i=1; i<=n; ++i)
	{
		ans += log10(i);
	}
	printf ("%d\n", (int)ans + 1);
}
/****************************************************************************/

/****************************************************************************/
计算几何

1. 基本函数

1.1 Point 定义
const double eps = 1e-8;
const double PI = acos(-1.0);
int sgn(double x)
{
	if(fabs(x) < eps)return 0;
	if(x < 0)return -1;
	else return 1;
}
struct Point
{
	double x,y;
	Point(){}
Point(double _x,double _y)
{
	x = _x;y = _y;
}
Point operator -(const Point &b)const
{
	return Point(x - b.x,y - b.y);
}
//叉积
double operator ^(const Point &b)const
{
	return x*b.y - y*b.x;
}
//点积
double operator *(const Point &b)const
{
	return x*b.x + y*b.y;
}
//绕原点旋转角度B（弧度值），后x,y的变化
void transXY(double B)
{
	double tx = x,ty = y;
	x = tx*cos(B) - ty*sin(B);
	y = tx*sin(B) + ty*cos(B);
}
};
1.2 Line 定义
struct Line
{
	Point s,e;
	Line(){}
Line(Point _s,Point _e)
{
	s = _s;e = _e;
}
//两直线相交求交点
//第一个值为0表示直线重合，为1表示平行，为0表示相交,为2是相交
//只有第一个值为2时，交点才有意义
pair<int,Point> operator &(const Line &b)const
{
	Point res = s;
	if(sgn((s-e)^(b.s-b.e)) == 0)
	{
	if(sgn((s-b.e)^(b.s-b.e)) == 0)
	return make_pair(0,res);//重合
	else return make_pair(1,res);//平行
	}
	double t = ((s-b.s)^(b.s-b.e))/((s-e)^(b.s-b.e));
	res.x += (e.x-s.x)*t;
	res.y += (e.y-s.y)*t;
	return make_pair(2,res);
}
};
1.3 两点间距离
//*两点间距离
double dist(Point a,Point b)
{
	return sqrt((a-b)*(a-b));
}
1.4 判断：线段相交
//*判断线段相交
bool inter(Line l1,Line l2)
{
	return
	max(l1.s.x,l1.e.x) >= min(l2.s.x,l2.e.x) &&
	max(l2.s.x,l2.e.x) >= min(l1.s.x,l1.e.x) &&
	max(l1.s.y,l1.e.y) >= min(l2.s.y,l2.e.y) &&
	max(l2.s.y,l2.e.y) >= min(l1.s.y,l1.e.y) &&
	sgn((l2.s-l1.e)^(l1.s-l1.e))*sgn((l2.e-l1.e)^(l1.s-l1.e)) <= 0 &&
	sgn((l1.s-l2.e)^(l2.s-l2.e))*sgn((l1.e-l2.e)^(l2.s-l2.e)) <= 0;
}
1.5 判断：直线和线段相交
//判断直线和线段相交
bool Seg_inter_line(Line l1,Line l2) //判断直线l1和线段l2是否相交
{
	return sgn((l2.s-l1.e)^(l1.s-l1.e))*sgn((l2.e-l1.e)^(l1.s-l1.e)) <= 0;
}
1.6 点到直线距离
//点到直线距离
//返回为result,是点到直线最近的点
Point PointToLine(Point P,Line L)
{
	Point result;
	double t = ((P-L.s)*(L.e-L.s))/((L.e-L.s)*(L.e-L.s));
	result.x = L.s.x + (L.e.x-L.s.x)*t;
	result.y = L.s.y + (L.e.y-L.s.y)*t;
	return result;
}
1.7 点到线段距离
//点到线段的距离
//返回点到线段最近的点
Point NearestPointToLineSeg(Point P,Line L)
{
	Point result;
	double t = ((P-L.s)*(L.e-L.s))/((L.e-L.s)*(L.e-L.s));
	if(t >= 0 && t <= 1)
	{
		result.x = L.s.x + (L.e.x - L.s.x)*t;
		result.y = L.s.y + (L.e.y - L.s.y)*t;
	}
	else
	{
		if(dist(P,L.s) < dist(P,L.e))
			result = L.s;
		else result = L.e;
	}
	return result;
}
1.8 计算多边形面积
//计算多边形面积
//点的编号从0~n-1
double CalcArea(Point p[],int n)
{
	double res = 0;
	for(int i = 0;i < n;i++)
	res += (p[i]^p[(i+1)%n])/2;
	return fabs(res);
}
1.9 判断点在线段上
//*判断点在线段上
bool OnSeg(Point P,Line L)
{
	return
	sgn((L.s-P)^(L.e-P)) == 0 &&
	sgn((P.x - L.s.x) * (P.x - L.e.x)) <= 0 &&
	sgn((P.y - L.s.y) * (P.y - L.e.y)) <= 0;
}
1.10 判断点在凸多边形内
//*判断点在凸多边形内
//点形成一个凸包，而且按逆时针排序（如果是顺时针把里面的<0改为>0）
//点的编号:0~n-1
//返回值：
//-1:点在凸多边形外
//0:点在凸多边形边界上
//1:点在凸多边形内
int inConvexPoly(Point a,Point p[],int n)
{
	for(int i = 0;i < n;i++)
	{
		if(sgn((p[i]-a)^(p[(i+1)%n]-a)) < 0)return -1;
		else if(OnSeg(a,Line(p[i],p[(i+1)%n])))return 0;
	}
	return 1;
}
1.11 判断点在任意多边形内
//*判断点在任意多边形内
//射线法，poly[]的顶点数要大于等于3,点的编号0~n-1
//返回值
//-1:点在凸多边形外
//0:点在凸多边形边界上
//1:点在凸多边形内
int inPoly(Point p,Point poly[],int n)
{
	int cnt;
	Line ray,side;
	cnt = 0;
	ray.s = p;
	ray.e.y = p.y;
	ray.e.x = -100000000000.0;//-INF,注意取值防止越界
	for(int i = 0;i < n;i++)
	{
		side.s = poly[i];
		side.e = poly[(i+1)%n];
		if(OnSeg(p,side))return 0;
		//如果平行轴则不考虑
		if(sgn(side.s.y - side.e.y) == 0)
			continue;
		if(OnSeg(side.s,ray))
		{
			if(sgn(side.s.y - side.e.y) > 0)cnt++;
		}
		else if(OnSeg(side.e,ray))
		{
			if(sgn(side.e.y - side.s.y) > 0)cnt++;
		}
		else if(inter(ray,side))
			cnt++;
	}
	if(cnt % 2 == 1)return 1;
	else return -1;
}
1.12 判断凸多边形
//判断凸多边形
//允许共线边
//点可以是顺时针给出也可以是逆时针给出
//点的编号1~n-1
bool isconvex(Point poly[],int n)
{
	bool s[3];
	memset(s,false,sizeof(s));
	for(int i = 0;i < n;i++)
	{
		s[sgn( (poly[(i+1)%n]-poly[i])^(poly[(i+2)%n]-poly[i]) )+1] = true;
		if(s[0] && s[2])return false;
	}
	return true;
}

2. 凸包
/*
* 求凸包，Graham算法
* 点的编号0~n-1
* 返回凸包结果Stack[0~top-1]为凸包的编号
*/
const int MAXN = 1010;
Point list[MAXN];
int Stack[MAXN],top;
//相对于list[0]的极角排序
bool _cmp(Point p1,Point p2)
{
	double tmp = (p1-list[0])^(p2-list[0]);
	if(sgn(tmp) > 0)return true;
	else if(sgn(tmp) == 0 && sgn(dist(p1,list[0]) - dist(p2,list[0])) <= 0)
	return true;
	else return false;
}
void Graham(int n)
{
	Point p0;
	int k = 0;
	p0 = list[0];
	//找最下边的一个点
	for(int i = 1;i < n;i++)
	{
		if( (p0.y > list[i].y) || (p0.y == list[i].y && p0.x > list[i].x) )
		{
			p0 = list[i];
			k = i;
		}
	}
	swap(list[k],list[0]);
	sort(list+1,list+n,_cmp);
	if(n == 1)
	{
		top = 1;
		Stack[0] = 0;
		return;
	}
	if(n == 2)
	{
		top = 2;
		Stack[0] = 0;
		Stack[1] = 1;
		return ;
	}
	Stack[0] = 0;
	Stack[1] = 1;
	top = 2;
	for(int i = 2;i < n;i++)
	{
		while(top > 1 && sgn((list[Stack[top-1]]-list[Stack[top-2]])
			^(list[i]-list[Stack[top-2]])) <=0)

		top--;
		Stack[top++] = i;
	}
}
/****************************************************************************/

/****************************************************************************/
黑科技

1. 快速读入输出(读入输出外挂)！
/*
	原理：将数字按照字符输出，加快速度
	适用场合：对象为整数(包括负数)，当读入大数据时才有效果
*/
inline int read(void)
{
    int x = 0, f = 1;	char ch = getchar ();
    while (ch < '0' || ch > '9')	{if (ch == '-')	f = -1;	ch = getchar ();}
    while (ch >= '0' && ch <= '9')	{x = x * 10 + ch - '0';	ch = getchar ();}
    return x * f;
}

inline void print(int x)
{
    if (x < 0)    {putchar ('-');	x = -x;}
    if (x > 9)    print (x / 10);
    putchar (x % 10 + '0');
}

2. 解决爆栈，手动加栈！
#pragma comment(linker, "/STACK:1024000000,1024000000")

3. #include <bits/stdc++.h>
/****************************************************************************/
