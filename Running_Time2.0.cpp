/****************************************************************************/
									目录
高精度
	1 大数模板-JayYe
	2 大数介乘求位数

排序
	1 快速排序
	2 归并排序

数学
	1 欧几里德算法
	2 扩展欧几里德算法
	3 乘方取模 && 乘法取模
	4 中国剩余定理
	5 逆元
	6 素数
	7 欧拉函数
	8 约瑟夫环公式
	9 第K个全排列
	10 矩阵快速幂
	
DP
	1 最长上升子序列
	2 最长公共子序列
	3 最长公共上升子序列
	4 最长回文串
	5 最大子序列和
	6 树形DP

字符串处理
	1 KMP
	2 Aho-Corasick
	3 Hash
	4 Manacher

数据结构
	1 线段树
	2 树状数组
	3 并查集
    4 莫队算法

图论
	1 最短路
		1.1 Dijkstra
		1.2 SPFA
		1.3 Floyd_Warshall
	2 拓扑排序
	3 最小生成树
		3.1 Kruskal
		3.2 Prim
	4 二分图最大匹配
	5 二分图最大权匹配
	6 LCA

表达式求值

黑科技
	1 快速读入输出(读入输出外挂)！
	2 解决爆栈，手动加栈！
	3 #include <bits/stdc++.h>
/****************************************************************************/

/****************************************************************************/
高精度

大数模板-JayYe
/*
    大数模版	-copy from JayYe
*/
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

const int numlen = 2005; // 需要的位数
const int numbit = 4;   // 数组一位表示的整数
const int addbit = 10000;//进位数
const int maxn = numlen/numbit + 10;   // 数组要开的位数

int Max(int a, int b) { return a>b?a:b; }
struct bign {
    int len, s[numlen];
    bign() {
        memset(s, 0, sizeof(s));
        len = 1;
    }
    bign(int num) { *this = num; }
    bign(const char *num) { *this = num; }
    bign operator = (const int num) {
        char s[numlen];
        sprintf(s, "%d", num);
        *this = s;
        return *this;
    }
    bign operator = (const char *num) {
        int clen = strlen(num);
        while(clen > 1 && num[0] == '0') num++, clen--;
        len = 0;
        for(int i = clen-1;i >= 0;i -= numbit) {
            int top = min(numbit, i+1), mul = 1;
            s[len] = 0;
            for(int j = 0;j < top; j++) {
                s[len] += (num[i-j]-'0')*mul;
                mul *= 10;
            }
            len++;
        }
        deal();
        return *this;
    }

    void deal() {
        while(len > 1 && !s[len-1]) len--;
    }

    bign operator + (const bign &a) const {
        bign ret;
        ret.len = 0;
        int top = Max(len, a.len) , add = 0;
        for(int i = 0;add || i < top; i++) {
            int now = add;
            if(i < len) now += s[i];
            if(i < a.len)   now += a.s[i];
            ret.s[ret.len++] = now%addbit;
            add = now/addbit;
        }
        return ret;
    }
    bign operator - (const bign &a) const {
        bign ret;
        ret.len = 0;
        int cal = 0;
        for(int i = 0;i < len; i++) {
            int now = s[i] - cal;
            if(i < a.len)   now -= a.s[i];
            if(now >= 0)    cal = 0;
            else {
                cal = 1; now += addbit;
            }
            ret.s[ret.len++] = now;
        }
        ret.deal();
        return ret;
    }
    bign operator * (const bign &a) const {
        bign ret;
        ret.len = len + a.len;
        for(int i = 0;i < len; i++) {
            int pre = 0;
            for(int j = 0;j < a.len; j++) {
                int now = s[i]*a.s[j] + pre;
                pre = 0;
                ret.s[i+j] += now;
                if(ret.s[i+j] >= addbit) {
                    pre = ret.s[i+j]/addbit;
                    ret.s[i+j] -= pre*addbit;
                }
            }
            if(pre) ret.s[i+a.len] = pre;
        }
        ret.deal();
        return ret;
    }
    //乘以小整数，直接乘快点  ***********注意计算过程可能会爆int
    bign operator * (const int num) {
        bign ret;
        ret.len = 0;
        int bb = 0;
        for(int i = 0;i < len; i++) {
            int now = bb + s[i]*num;
            ret.s[ret.len++] = now%addbit;
            bb = now/addbit;
        }
        while(bb) {
            ret.s[ret.len++] = bb % addbit;
            bb /= addbit;
        }
        ret.deal();
        return ret;
    }
    // 除以一个小整数     ***********注意计算过程可能会爆int
    bign operator / (const int a) const {
        bign ret;
        ret.len = len;
        int pre = 0;
        for(int i = len-1;i >= 0; i--) {
            ret.s[i] = (s[i] + pre*addbit)/a;
            pre = s[i] + pre*addbit - a*ret.s[i];
        }
        ret.deal();
        return ret;
    }
    bign operator % (const int a) const {
        bign b = *this / a;
        return *this - b*a;
    }

    bign operator += (const bign &a) { *this = *this + a; return *this; }
    bign operator -= (const bign &a) { *this = *this - a; return *this; }
    bign operator *= (const bign &a) { *this = *this * a; return *this; }
    bign operator /= (const int a) { *this = *this / a; return *this; }
    bign operator %= (const int a) { *this = *this % a; return *this; }

    bool operator < (const bign &a) const {
        if(len != a.len)    return len < a.len;
        for(int i = len-1;i >= 0; i--) if(s[i] != a.s[i])
            return s[i] < a.s[i];
        return false;
    }
    bool operator > (const bign &a) const  { return a < *this; }
    bool operator <= (const bign &a) const { return !(*this > a); }
    bool operator >= (const bign &a) const { return !(*this < a); }
    bool operator == (const bign &a) const { return !(*this > a || *this < a); }
    bool operator != (const bign &a) const { return *this > a || *this < a; }
};
istream& operator >> (istream &in, bign &x) {
	string s;
	in >> s;
	x = s.c_str();
	return in;
}
ostream& operator << (ostream &out, const bign &x) {
    printf("%d", x.s[x.len-1]);
    for(int i = x.len-2;i >= 0; i--)    printf("%04d", x.s[i]);
    return out;
}

bign f[110];	//卡特兰数
int main(void)	{
	f[0] = 1;
	for (int i=1; i<=100; ++i)	{
		f[i] = f[i-1] * (4 * i - 2) / (i + 1);
	}
	int n;
	while (scanf ("%d", &n) == 1)	{
		if (n == -1)	break;
		cout << f[n] << endl;
	}

	return 0;
}

2 大数介乘求位数
/*
	log10(n!) = log10(1*2*..*n) = log10(1) + log10(2) + ...+log10(n)
	解释：123456=1.23456*10^5;
		  log10(123456)=5.09151;
		  log10(1.23456*10^5)=log10(1.23456)+log10(10^5)=0.09151+5;
		  故int(log10(n))+1 就是n的位数
*/
void work(int n)	{
	double ans = 0;
	for (int i=1; i<=n; ++i)	{
		ans += log10(i);
	}
	printf ("%d\n", (int)ans + 1);
}
/****************************************************************************/

/****************************************************************************/
排序

1 快速排序
/*
	quick_sort：快速排序，随机标杆，O(nlogn)
*/
void quick_sort(int *a, int l, int r)   {
    if (l < r)  {
        swap (a[l+rand()%(r-l)], a[l]);
        int i = l, j = r, x = a[l];
        while (i < j)   {
            while (i < j && a[j] >= x)  j--;
            if (i < j)  a[i++] = a[j];
            while (i < j && a[i] < x)   i++;
            if (i < j)  a[j--] = a[i];
        }
        a[i] = x;
        quick_sort (a, l, i-1);
        quick_sort (a, i+1, r);
    }
}
/*
	基于快速排序思想，返回第K小的数字
*/
int Kth_small(int *a, int l, int r, int k)	{
	if (l == r)	return a[l];
	int x = a[l+rand ()%(r-l+1)], i = l - 1, j = r + 1;
	while (i < j)	{
		while (a[++i] < x);	while (a[--j] > x);
		if (i < j)	swap (a[i], a[j]);
	}
	if (j == r)	j--;
	i = j - l + 1;
	if (k <= i)	return Kth_small (a, l, j, k);
	else	return Kth_small (a, j+1, r, k - i);
}

2 归并排序
/*
	merge_sort：归并排序，O(nlogn)
*/
void merge(int *a, int p, int q, int r)   {
    int n1 = q - p + 1, n2 = r - q;
    for (int i=1; i<=n1; ++i)   L[i] = a[p+i-1];
    for (int i=1; i<=n2; ++i)   R[i] = a[q+i];
    L[n1+1] = R[n2+1] = INF;
    for (int i=1, j=1, k=p; k<=r; ++k)  {	//求逆序数：cnt += n1 - i + 1;
        a[k] = (L[i] <= R[j]) ? L[i++] : R[j++];
    }
}
void merge_sort(int *a, int p, int r)   {
    if (p < r)  {
        int q = (p + r) >> 1;
        merge_Sort (a, p, q);
        merge_Sort (a, q+1, r);
        merge (a, p, q, r);
    }
}
/****************************************************************************/

/****************************************************************************/
数学

1 欧几里德算法
/*
	辗转相除法(欧几里德算法)
*/
int GCD(int a, int b)   {
    return b ? GCD (b, a % b) : a;
}
/*
	非递归版本
*/
int GCD2(int a, int b)  {
    while (b)   {
        int c = b;
        b = a % b;
        a = c;
    }
    return a;
}
/*
	快速GCD
*/
int quick_GCD(int a, int b)  {
    if (!a) return b;
    if (!b) return a;
    if (! (a & 1) && ! (b & 1)) return quick_GCD (a >> 1, b >> 1) << 1;
    else if (! (b & 1)) return quick_GCD (a, b >> 1);
    else if (! (a & 1)) return quick_GCD (a >> 1, b);
    else    return quick_GCD (abs (a - b), min (a, b));
}

2 扩展欧几里德算法
/*
	扩展欧几里得算法， 求x, y 使得 GCD(a, b) = a * x + b * y;
	返回d = GCD (a,b)，和对应于等式 ax + by = d 中的x，y
*/
int ex_GCD(int a, int b, int &x, int &y)    {
    if (!a && !b)   return -1;					//无最大公约数
    if (!b) {
        x = 1;  y = 0;  return a;
    }
    int d = ex_GCD (b, a % b, y, x);
    y -= a / b * x;
    return d;
}
/*
	模线性方程 a * x = b (% n)	By 吉林模板
*/
void modeq(int a, int b, int n) {	// ! n > 0
    int e, i, d, x, y;
    d = extgcd(a, n, x, y);
    if (b % d > 0) printf("No answer!\n");
    else	{
        e = (x * (b / d)) % n;
        for (i = 0; i < d; i++) // !!! here x maybe < 0
            printf("%d-th ans: %d\n", i+1, (e+i*(n/d))%n);
    }
}

3 乘方取模 && 乘法取模
/*
	乘方取模(快速幂)，复杂度O (logn)
*/
int pow_mod(int x, int n, int p)    {       //x ^ n % p
    int ret = 1;
    while (n)   {
        if (n & 1)  {
            ret = 1ll * ret * x % p;
        }
        x = 1ll * x * x % p;
        n >>= 1;
    }
    return ret;
}
/*
	乘法取模，复杂度O (logb)
*/
int multi_mod(int a, int b, int p)  {       //a * b % p
    int ret = 0;
    a = (a % p + p) % p;
    b = (b % p + p) % p;
    while (b)   {
        if (b & 1)  {
            ret += a;
            if (ret >= p)   ret -= p;
        }
        b >>= 1;
        a <<= 1;
        if (a >= p) a -= p;
    }
    return ret;
}

4 中国剩余定理
/*
	求解满足：x = 2 (%3), x = 3 (%5), x = 2 (%7)的最小x
	解为：x = 2 * 35a + 3 * 21b + 2 * 15c (%105)
	k是方程组的个数，b对应2，3，2，m对应3，5，7
*/
int China(int k, int *b, int *m) {
    ll M = 1, x, y, ret = 0;
    for (int i=1; i<=k; ++i)    M *= m[i];
    for (int i=1; i<=k; ++i)    {
        ll w = M / m[i];
        ex_GCD (w, m[i], x, y);
        ret += multi_mod (multi_mod (x, w, M), b[i], M);		//防爆long long
    }
    return (ret + M) % M;
}

5 逆元
/*
	求法1：预处理，数组保存，复杂度 O (n)
*/
void Inv(int n, int p)    {                        //p是质数
    inv[1] = 1;
    for (int i=2; i<=n; ++i)    {
        inv[i] = 1ll * (p - p / i) * inv[p%i] % p;
    }
}
/*
	求法2：复杂度 O (logn)
*/
int Inv(int a, int p)   {                           //p是质数
    return pow_mod (a, p - 2, p);
}
/*
	求法3：扩展欧几里德求解
*/
int Inv(int a, int p)  {	                        //p是质数
    int x, y;
    int d = ex_GCD (a, p, x, y);
    if (d == 1) return (x % p + p) % p;
    else	return -1;
}

6 素数
/*
    欧拉筛法：返回n以内素数的个数， 复杂度O (n)
*/
int seive(int n)    {        //Euler (欧拉筛法)
    int p = 0;
    memset (is_prime, true, sizeof (is_prime));
    is_prime[0] = is_prime[1] = false;
    for (int i=2; i<=n; ++i)    {
        if (is_prime[i])    prime[++p] = i;
        for (int j=1; j<=p && i*prime[j]<=n; ++j)    {
            is_prime[i*prime[j]] = false;
            if (i % prime[j] == 0)    break;
        }
    }
    return p;
}

/*
    素性测试，在小范围(1e5)内判素数个数以及单个数判素数有奇效
    不适用于大范围判素数
*/
bool is_prime(int x)    {
    if (x == 2 || x == 3)   return true;
    if (x % 6 != 1 && x % 6 != 5)   return false;
    for (int i=5; i*i<=x; i+=6) {
        if (x % i == 0 || x % (i + 2) == 0) return false;
    }
    return true;
}
/*
	素性测试，is_prime[x] = true / false，不能得到prime[i]
*/
void seive(int n)   {
    memset (is_prime, true, sizeof (is_prime));
    is_prime[0] = is_prime[1] = false;
    int m = sqrt (n);
    for (int i=2; i<=m; ++i)    {
        if (is_prime[i])    {
            for (int j=i*i; j<=n; j+=i)   {
                is_prime[j] = false;
            }
        }
    }
}
/*
	素性测试，Miller_Rabin 随机算法
    以a为基,n-1=x*2^t，a^(n-1) = 1(mod n)  验证n是不是合数
    一定是合数返回true, 不一定返回false
*/
bool Miller_Rabin(ll n) {
    if (n == 2) return true;
    if (n < 2 || ! (n & 1))  return false;           //偶数或1
    ll x = n - 1;   int t = 0;
    while (! (x & 1))   {
        x >>= 1;  t++;
    }
    for (int i=1; i<=S; ++i) {						//srand (time (NULL));
        ll a = rand () % (n - 1) + 1;               //需要cstdlib，ctime头文件
        if (check (a, n, x, t)) return false;       //合数
    }
    return true;
}
bool check(ll a, ll n, ll x, int t) {
    ll ret = pow_mod (a, x, n);
    ll last = ret;
    for (int i=1; i<=t; ++i) {
        ret = multi_mod (ret, ret, n);
        if (ret == 1 && last != 1 && last != n - 1) return true;    //合数
        last = ret;
    }
    if (ret != 1)   return true;
    return false;
}

/*
    整数分解，先打个素数表优化试除法
*/
vector<int> factorize(int n)  {
    int p = seive (100000);
    vector<int> ret;
    for (int i=1; i<=p && prime[i]<=n/prime[i]; ++i)    {
        while (n % prime[i] == 0)  {
            n /= prime[i];
            ret.push_back (prime[i]);
        }
    }
    if (n != 1) ret.push_back (n);
    return ret;
}
/*
    大整数分解，Pollard_rho 随机算法
    factorize ()保存质因数在vector
*/
ll Pollard_rho(ll x, ll c)  {
    ll i = 1, k = 2;
    ll a = rand () % x;
    ll b = a;
    while (1)   {
        i++;
        a = (multi_mod (a, a, x) + c) % x;
        ll d = GCD (b - a, x);
        if (d != 1 && d != x)   return d;
        if (b == a) return x;
        if (i == k) b = a, k += k;
    }
}
void factorize(ll n, vector<ll> &ret) {           //ret保存质因数，无序
    if (Miller_Rabin (n))   {                       //素数
        ret.push_back (n);  return ;
    }
    ll p = n;
    while (p >= n)   p = Pollard_rho (p, rand () % (n - 1) + 1);
    factorize (p, ret);
    factorize (n / p, ret);
}

/*
	约数枚举，复杂度O (sqrt(n))
	By TiaoZhan
*/
vector<int> divisor(int n)	{
    vector<int> res;
    for (int i=1; i*i<=n; ++i)	{
		if (n % i == 0)	{
			res.push_back (i);
			if (n / i != i)	res.push_back (n / i);
		}
	}
	return res;
}

7 欧拉函数
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
unsigned euler(unsigned x)	{
    //  就是公式
    unsigned i, res=x;
    for (i = 2; i < (int)sqrt(x * 1.0) + 1; i++)
        if (x%i==0)	{
            res = res / i * (i - 1);
            while (x % i == 0) x /= i; //  保证i 一定是素数
        }
    if (x > 1) res = res / x * (x - 1);
    return res;
}

8 约瑟夫环公式
/*
	经典的约瑟夫(Joseph)水题
	(转载)约瑟夫环的递推公式：
    1. f[1]=0;　f[i]=(f[i-1]+m)%i; (i>1)	(0~n-1)
    2. f[1]=1;　f[i]=(f[i-1]+m)%i  (i>1);   if(f[i]==0) f[i]=i;		(1~n)
    3. P(1, m, k)=1 (i = 1);   P(i, m, k)=[P(i - 1, m, k ) + m - 1] % i + 1； 
		(i > 1, 此处先减1是为了让模i的值不为0)	(1~n)
*/

9 第K个全排列
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
long long fact(int x)	{
    long long res = 1;
    for (int i=1; i<=x; ++i)    res *= i;
    return res;
}
long long f(int step)	{
    long long res = fact (step);
    for (int i=0; i<26; ++i)    if (cnt[i])    res /= fact (cnt[i]);
    return res;    
}
void DFS(int step, long long k)	{
    if (step == len)	{
        for (int i=0; i<len; ++i)    printf ("%c", 'A' + pos[i]);
        puts ("");    return ;
    }

    for (int i=0; i<26; ++i)	{
        if (cnt[i] == 0)    continue;
        cnt[i]--;
        long long tmp = f (len - step - 1);
        if (tmp < k)    {k -= tmp;    cnt[i]++;}
        else	{
            pos[step] = i;
            DFS (step+1, k);
            return ;
        }
    }
}
10 矩阵快速幂
/*
	矩阵快速幂：求第n项的Fibonacci数，转置矩阵都给出，套个模板就可以了。效率很高啊
	Mat x;
    x.m[0][0] = x.m[0][1] = x.m[1][0] = 1; x.m[1][1] = 0;
    printf ("%d\n", matrix_pow_mod (x, n));
*/
struct Mat   {
    int m[2][2];
};
Mat multi_mod(Mat a, Mat b)   {
    Mat ret;    memset (ret.m, 0, sizeof (ret.m));
    for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
            for (int k=0; k<2; ++k) {
                ret.m[i][j] = (ret.m[i][j] + a.m[i][k] * b.m[k][j]) % MOD;
            }
        }
    }
    return ret;
}
int matrix_pow_mod(Mat x, int n)   {
    Mat ret; ret.m[0][0] = ret.m[1][1] = 1;  ret.m[0][1] = ret.m[1][0] = 0;
    while (n)   {
        if (n & 1)  ret = multi_mod (ret, x);
        x = multi_mod (x, x);
        n >>= 1;
    }
    return ret.m[0][1];
}
/*
    POJ_3735 题意：k次操作，g：i猫+1, e：i猫eat，s：swap
    矩阵快速幂：写个转置矩阵，将k次操作写在第0行，
    定义A = {1，0, 0, 0...}除了第一个外其他是猫的初始值	
*/
typedef long long ll;
const int MAXN = 1e2 + 10;
const int INF = 0x3f3f3f3f;
struct Mat  {
    ll m[N][N];
    Mat ()  {
        memset (m, 0, sizeof (m));
    }
    void init(void) {
        for (int i=0; i<N; ++i)  m[i][i] = 1;
    }
};
int n, m, k;
Mat operator * (Mat &a, Mat &b) {
    Mat ret;
    for (int k=0; k<=n; ++k)    {
        for (int i=0; i<=n; ++i)    {
            if (a.m[i][k])  {
                for (int j=0; j<=n; ++j)    {
                    ret.m[i][j] += a.m[i][k] * b.m[k][j];
                }
            }
        }
    }
    return ret;
}
Mat operator ^ (Mat x, int n)   {
    Mat ret;    ret.init ();
    while (n)   {
        if (n & 1)  ret = ret * x;
        x = x * x;
        n >>= 1;
    }
    return ret;
}
int main(void)  {       //POJ 3735 Training little cats
    while (scanf ("%d%d%d", &n, &m, &k) == 3)   {
        if (!n && !m && !k) break;
        Mat A, T;   A.m[0][0] = 1;  T.init ();
        char op[5]; int p, q;
        while (k--) {
            scanf ("%s", op);
            if (op[0] == 'g')    {
                scanf ("%d", &p);   T.m[0][p]++;
            }
            else if (op[0] == 'e')  {
                scanf ("%d", &p);
                for (int i=0; i<=n; ++i)    T.m[i][p] = 0;
            }
            else    {
                scanf ("%d%d", &p, &q);
                for (int i=0; i<=n; ++i)    swap (T.m[i][p], T.m[i][q]);
            }
        }
        Mat ans = A * (T ^ m);
        for (int i=1; i<=n; ++i)    printf ("%I64d%c", ans.m[0][i], (i==n) ? '\n' : ' ');
    }

    return 0;
}
/****************************************************************************/

/****************************************************************************/
DP

1 最长上升子序列
/*
    LIS(Longest Increasing Subsequence)  二分查找优化, O(nlogn)
    设当前最长递增子序列为len，考虑元素a[i]; 若d[len]<a[i]，则len++，并使d[len]=a[i];
    否则,在d[1~len]中二分查找第一个大于等于a[i]的位置j，使d[j]=a[i]。附上打印路径代码(准确性未知)
*/
void LIS(void)  {
    int len = 1;    d[1] = a[1];    fa[1] = -1;
    for (int i=2; i<=n; ++i) {
        if (d[len] < a[i])   {
            d[++len] = a[i];
            pos[len] = i;   fa[i] = pos[len-1];
        }
        else    {
            int j = lower_bound (d+1, d+1+len, a[i]) - d;
            d[j] = a[i];
            pos[j] = i; fa[i] = (j == 1) ? -1 : pos[j-1];
        }
    }
    printf ("%d\n", len);
    vector<int> res;  int i;
    for (i=pos[len]; ~fa[i]; i=fa[i])   res.push_back (a[i]);
    res.push_back (a[i]);
    for (int i=res.size ()-1; i>=0; --i) printf ("%d%c", res[i], i == 0 ? '\n' : ' ');
}

2 最长公共子序列
/*
    LCS(Longest Common Subsequence)：
        状态转移方程：dp[i][j] = dp[i-1][j-1] + 1; (s[i-1] == t[i-1])
                    dp[i][j] = max (dp[i][j-1], dp[i-1][j]);(s[i-1] != t[i-1])
        可滚动数组优化。附带有print输出路径函数。
*/
void LCS(void)  {
    memset (dp, 0, sizeof (dp));
    memset (fa, 0, sizeof (fa));
    for (int i=0; i<=lens; ++i)    fa[i][0] = -1;
    for (int i=0; i<=lent; ++i)    fa[0][i] = 1;
 
    for (int i=1; i<=lens; ++i) {
        for (int j=1; j<=lent; ++j) {
            if (s[i-1] == t[j-1])   {
                dp[i][j] = dp[i-1][j-1] + 1;
                fa[i][j] = 0;
            }
            else if (dp[i-1][j] >= dp[i][j-1])  {
                dp[i][j] = dp[i-1][j];
                fa[i][j] = -1;
            }
            else    {
                dp[i][j] = dp[i][j-1];
                fa[i][j] = 1;
            }
        }
    }
 
    printf ("%d\n", dp[lens][lent]);
    print (lens, lent); puts ("");
}
void print(int x, int y)    {
    if (!x && !y)   return ;
    if (fa[x][y] == 0)  {
        print (x-1, y-1);    printf ("%c", s[x-1]);
    }
    else if (fa[x][y] == -1)    {
        print (x-1, y);    printf ("%c", s[x-1]);
    }
    else    {
        print (x, y-1);    printf ("%c", t[y-1]);
    }
}

3 最长公共上升子序列
/*
    LCIS(Longest Common Increasing Subsequence) 最长公共上升子序列
    状态转移方程：a[i] != b[j]: dp[i][j] = dp[i-1][j];
                  a[i] == b[j]: dp[j]=max(dp[j],dp[k]); (1<=k<j&&b[k]<b[j])　　
    打印路径时按照b[i]来输出
*/
void LCIS(void) {
    memset (dp, 0, sizeof (dp));
    memset (fx, 0, sizeof (fx));
    memset (fy, 0, sizeof (fy));
    int sx = 0, sy = 0;
    int ret = 0, k = 0;
    for (int i=1; i<=n; ++i) {
        k = 0;
        for (int j=1; j<=m; ++j) {
            dp[i][j] = dp[i-1][j];               //以a[]为主循环，每个a[i]，去找每个b[j]
            fx[i][j] = i - 1;   fy[i][j] = j;
            if (a[i] == b[j] && dp[i][j] < dp[i][k] + 1) {		  //满足LCS
                dp[i][j] = dp[i][k] + 1;        //在1~j-1找到b[k]<a[i]，满足LIS，在b[k]上更新dp
                fx[i][j] = i;   fy[i][j] = k;
            }
            else if (a[i] > b[j] && dp[i][j] > dp[i][k])  k = j;  //找到最优的k
            if (ret < dp[i][j])  {
                ret = dp[i][j];                 //更新所有dp中的最大值
                sx = i, sy = j;
            }
        }
    }
    printf ("%d\n", ret);
    fir = true;
    print (sx, sy, -1); puts ("");
}	
void print(int x, int y, int last)  {       //bool fir;
    if (x == 0 || y == 0)   return ;
    print (fx[x][y], fy[x][y], y);
    if (y != last)  {
        if (fir)    printf ("%d", b[y]), fir = false;
        else    printf (" %d", b[y]);
    }
}

4 最长回文串
/*
    LCS的思想，dp[i][j]表示i到j的最长回文串长度，状态转移方程：
        1. dp[j][j+i-1] = dp[j+1][j+i-2] + 2;   (str[j] == str[j+i-1])
        2. dp[j][j+i-1] = max (dp[j+1][j+i-1], dp[j][j+i-2]);   (str[j] != str[j+i-1])
    ans[1][len]是string类型，记录LPS字符(Longest Palidromic Subsequence)
*/
void LPS(void)  {
    int len = strlen (str + 1);
    memset (dp, 0, sizeof (dp));
    for (int i=1; i<=len; ++i)  dp[i][i] = 1;
    for (int i=1; i<=len; ++i)  ans[i][i] = str[i];
 
    for (int i=2; i<=len; ++i)   {              //区间长度
        for (int j=1; j+i-1<=len; ++j)  {       //[j, j+i-1]
            if (str[j] == str[j+i-1])   {
                if (i == 2) {
                    dp[j][j+i-1] = 2;
                    ans[j][j+i-1] = ans[j][j] + ans[j+i-1][j+i-1];  continue;
                }
                dp[j][j+i-1] = dp[j+1][j+i-2] + 2;
                ans[j][j+i-1] = str[j] + ans[j+1][j+i-2] + str[j+i-1];
            }
            else if (dp[j+1][j+i-1] > dp[j][j+i-2]) {
                dp[j][j+i-1] = dp[j+1][j+i-1];
                ans[j][j+i-1] = ans[j+1][j+i-1];
            }
            else if (dp[j][j+i-2] > dp[j+1][j+i-1]) {
                dp[j][j+i-1] = dp[j][j+i-2];
                ans[j][j+i-1] = ans[j][j+i-2];
            }
            else    {
                dp[j][j+i-1] = dp[j+1][j+i-1];
                ans[j][j+i-1] = min (ans[j+1][j+i-1], ans[j][j+i-2]);
            }
        }
    }
    int mlen = dp[1][len];
    for (int i=0; i<mlen; ++i)  {
        printf ("%c", ans[1][len][i]);
    }
    puts ("");
}

5 最大子序列和
/*
    MCS (Maximum Continuous Subsequence) 最大子序列和 O (n)
    1. DP   2. 前缀
    若有多个答案输出第一个，均给出区间端点
*/
void MCS(int n) {
    int l = 0, ll = 0, rr = 0;
    int sum = -INF, mx = -INF;
    for (int i=1; i<=n; ++i) {
        if (sum + a[i] < a[i])   {
            sum = a[i]; l = i;
        }
        else    sum += a[i];
        if (sum > mx)    {
            mx = sum;   ll = l, rr = i;
        }
    }
    printf ("%d %d %d\n", mx, ll, rr);
}
//O (n) //another
void MCS(int n) {
    int l = 0, ll = 0, rr = 0;
    int sum = 0, mx = -INF, mn = 0;
    for (int i=1; i<=n; ++i) {
        sum += a[i];
        if (sum - mn > mx)   {
            mx = sum - mn;  ll = l; rr = i;
        }
        if (sum < mn)    {
            mn = sum;   l = i;
        }
    }
    printf ("%d %d %d\n", mx, ll + 1, rr);
}

6 树形DP
/*
    URAL_1039 题意：上司在，员工不在，反之不一定。每一个人有一个权值，问权值和最大多少。
    树形DP：把上司和员工的关系看成根节点和子节点的关系，两者有状态转移方程：
			dp[rt][0] += max (dp[son][1], dp[son][0]);	//上司不去
			dp[rt][1] += dp[son][0];					//上司去，员工都不去
*/
void DFS(int rt)    {
	vis[rt] = true;
	dp[rt][0] = 0;
	for (int i=0; i<edge[rt].size (); ++i)  {
		int son = edge[rt][i];
		if (!vis[son])  {
			DFS (son);
			dp[rt][0] += max (dp[son][1], dp[son][0]);
			dp[rt][1] += dp[son][0];
		}
	}
}
/****************************************************************************/

/****************************************************************************/
字符串处理
1 KMP
/*
    T[]是文本串，P[]是模式串，lent，lenp是各自的长度
    返回第一个匹配完成在文本串上的位置，不成功返回-1
*/
void get_fail(char *P, int lenp)  {
    int i = 0, j = -1;  fail[0] = -1;
    while (i < lenp)    {
        if (j == -1 || P[j] == P[i])    {
            i++;    j++;    fail[i] = j;
        }
        else    j = fail[j];
    }
}
int KMP(char *T, char *P)   {
    int lent = strlen (T), lenp = strlen (P);
    get_fail (P, lenp);
    int i = 0, j = 0;
    while (i < lent)  {
        while (j != -1 && T[i] != P[j]) j = fail[j];
        i++; j++;
        if (j == lenp)  return (i - j + 1);
    	// if (j >= lenp) {						//求出现的个数
        	// ans++;  j = nex[j];				//可重复
    	// }									//j = 0重新匹配
    }
    return -1;
}

2 Aho-Corasick
/*
    HDOJ 2222 题意：每个文本串的出现次数
    AC自动机：入门题，注意关键词重复和出现过的关键词的次数不算到总数内   
*/
const int MAXN = 1e4 + 10;
const int MAXNODE = MAXN * 50 + 10;
const int SIZE = 26;
struct AC   {
    int ch[MAXNODE][SIZE], fail[MAXNODE], val[MAXNODE], sz;
    void init(void) {
        memset (ch[0], 0, sizeof (ch[0]));
        sz = 1; val[0] = 0;
    }
    int idx(char c) {
        return c - 'a';
    }
    void insert(char *P)    {
        int u = 0;
        for (int i=0; P[i]; ++i)   {
            int c = idx (P[i]);
            if (!ch[u][c])  {
                memset (ch[sz], 0, sizeof (ch[sz]));
                ch[u][c] = sz;  val[sz++] = 0;
            }
            u = ch[u][c];
        }
        val[u]++;
    }
    void get_fail(void) {
        queue<int> Q;   fail[0] = 0;
        for (int i=0; i<SIZE; ++i)  {
            int u = ch[0][i];
            if (u)  {
                fail[u] = 0;    Q.push (u);
            }
        }
        while (!Q.empty ()) {
            int u = Q.front (); Q.pop ();
            for (int i=0; i<SIZE; ++i)  {
                int &v = ch[u][i];
                if (!v) {
                    v = ch[fail[u]][i];  continue;
                }
                Q.push (v);
                fail[v] = ch[fail[u]][i];   //val[v] += val[fail[u]];
            }
        }
    }
    int query(char *T)  {
        int ret = 0;
        for (int u=0, i=0; T[i]; ++i)   {
            int c = idx (T[i]);
            u = ch[u][c];
            ret += val[u];
            int tmp = u;
            while (tmp) {
                ret += val[tmp];    val[tmp] = 0;
                tmp = fail[tmp];
            }
        }
        return ret;
    }
}ac;
char pat[55], tex[1000010];
int main(void)    {
    int T;  scanf ("%d", &T);
    while (T--) {
        int n;  scanf ("%d", &n);
        ac.init ();
        for (int i=1; i<=n; ++i)    {
            scanf ("%s", &pat);   ac.insert (pat);
        }
        ac.get_fail (); scanf ("%s", &tex);
        printf ("%d\n", ac.query (tex));
    }
    return 0;
}

3 Hash
/*
	Rabin_Karp：字符串hash，RK算法的原理：首先把模式串的哈希值算出来，
    在文本串里不断更新模式串的长度的哈希值，若相等，则找到了，否则整个模式串的长度的哈希值向右移动一位
	给你两个字符串A，B，请输出B字符串在A字符串中出现了几次。
*/
const ull KEY = 100000007;
int Rabin_Karp(void) {
    int lent = strlen (t);
    int lenp = strlen (p);
    ull h = 1;
    for (int i=0; i<lenp; ++i)    h *= KEY;
    ull th = 0, ph = 0; int ret = 0;
    for (int i=0; i<lenp; ++i)    th = th * KEY + t[i];
    for (int i=0; i<lenp; ++i)    ph = ph * KEY + p[i];
    for (int i=0; i+lenp<=lent; ++i)    {
        if (th == ph)   {
            ret++;
            for (int j=i; j<=i+lenp-1; ++j) {       //找到了一个模式串，不能再用，整个跳过去
                th = th * KEY + t[j+lenp] - t[j] * h;
            }
            i += lenp - 1;  continue;
        }
        if (i + lenp < lent) {
            th = th * KEY + t[i+lenp] - t[i] * h;
        }
    }
    return ret;
}

4 Manacher
/*
	最长回文串算法
*/
int Manacher(void)  {
    int len = strlen (s);
    for (int i=len; i>=0; --i)  {
        s[i*2+2] = s[i];    s[i*2+1] = '#';
    }
    s[0] = '$'; len = len * 2 + 2;
    int id = 0; p[0] = 1;
    for (int i=2; i<len; ++i)   {
        if (id + p[id] > i) p[i] = min (p[2*id-i], id + p[id] - i);
        else    p[i] = 1;
        while (s[i-p[i]] == s[i+p[i]])  p[i]++;
        if (id + p[id] < i + p[i])  id = i;
    }
    int ret = 0;
    for (int i=0; i<len-1; ++i)   {
        if (ret < p[i]) ret = p[i] - 1;
    }
    return ret;
}
/****************************************************************************/

/****************************************************************************/
数据结构

1 线段树
/*
	线段树：成段更新，延迟标记，即每次更新不更新到底，延迟到下一次更新or查询
	以下代码完成区间修改值以及区间查询和
*/
#define lson l, mid, rt << 1
#define rson mid + 1, r, rt << 1 | 1
struct ST   {
    int sum[N<<2], col[N<<2];
    void push_up(int rt)    {
        sum[rt] = sum[rt<<1] + sum[rt<<1|1];
    }
    void push_down(int rt, int len)   {
        if (col[rt])    {
            col[rt<<1] = col[rt<<1|1] = col[rt];
            sum[rt<<1] += col[rt] * (len - (len >> 1));
            sum[rt<<1|1] += col[rt] * (len >> 1);
            col[rt] = 0;
        }
    }
    void build(int l, int r, int rt)    {
        col[rt] = 0;
        if (l == r) {
            scanf ("%d", &sum[rt]);    return ;
        }
        int mid = (l + r) >> 1;
        build (lson); build (rson);
        push_up (rt);
    }
    void updata(int ql, int qr, int c, int l, int r, int rt)    {
        if (ql <= l && r <= qr) {
            sum[rt] += c * (r - l + 1); col[rt] += c;   return ;
        }
        push_down (rt, r - l + 1);
        int mid = (l + r) >> 1;
        if (ql <= mid)  updata (ql, qr, c, lson);
        if (qr > mid)   updata (ql, qr, c, rson);
        push_up (rt);
    }
    int query(int ql, int qr, int l, int r, int rt)    {
        if (ql <= l && r <= qr) return sum[rt];
        push_down (rt, r - l + 1);
        int mid = (l + r) >> 1, ret = 0;
        if (ql <= mid)  ret += query (ql, qr, lson);
        if (qr > mid)   ret += query (ql, qr, rson);
        return ret;
    }
}st;
/*
	区间合并：lsum[]统计从左端点起最长连续空房间数，rsum[]类似，sum[]统计区间最长连续的空房间数，
    它有三种情况：1.纯粹是左端点起的房间数；2.纯粹是右端点的房间数；3.当从左(右)房间起都连续时，
    加上另一个子节点，从左(右)房间起的数，sum[]再求最大值更新维护。
*/
void push_up(int rt, int len)    {
    lsum[rt] = lsum[rt<<1];
    rsum[rt] = rsum[rt<<1|1];
    if (lsum[rt] == len - (len>>1)) lsum[rt] += lsum[rt<<1|1];
    if (rsum[rt] == len>>1)   rsum[rt] += rsum[rt<<1];
    sum[rt] = max (rsum[rt<<1] + lsum[rt<<1|1], max (sum[rt<<1], sum[rt<<1|1]));
}


2 树状数组
/*
	以下代码实现功能与线段树同，sum (l, r) = query (r) - query (l-1)
*/
struct BIT  {
    int c[N], SZ;
    void init(int n) {
        memset (c, 0, sizeof (c));
        SZ = n;
    }
    void updata(int i, int x)   {
        while (i <= SZ)   {
            c[i] += x;  i += i & (-i);
        }
    }
    int query(int i)    {
        int ret = 0;
        while (i)   {
            ret += c[i];    i -= i & (-i);
        }
        return ret;
    }
}bit;


3 并查集
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

4 莫队算法
/*
    例题HDOJ_4358: 给你一棵树，树上的每个节点都有树值，给m个查询，问以每个点u为根的子树下有多少种权值恰好出现k次。
分析：首先要对权值离散化，然后要将树形转换为线形。然后按照右端点从小到大排序，
离线操作：将每一个深度的权值分组到相同权值的cnt中，当sz == k时，用树状数组更新+1，表示在该深度已经存在k个相同的权值，
如果>k，之前k个-2(-1是恢复原样，再-1是为下次做准备？)，然后一个点的子树的答案就是 sum (r) - sum (l-1)。
*/
const int N = 1e5 + 10;
const int INF = 0x3f3f3f3f;
const int MOD = 1e9 + 7;
struct Edge {
    int v, nex;
}edge[N*2];
struct Data {
    int l, r, id, b;
    Data () {}
    Data (int l, int r, int id, int b) : l (l), r (r), id (id), b (b) {}
    bool operator < (const Data &x) const {
        if (b == x.b)   return r < x.r;
        else    return b < x.b;
    }
}data[N];
int head[N], dfn[N], low[N], w[N], p[N], val[N], ans[N];
int cnt[N];
int n, k, m, e, dep, sum;
  
void init(void) {
    memset (head, -1, sizeof (head));
    e = 0;  dep = 0;
}
  
void add_edge(int u, int v) {
    edge[e].v = v;  edge[e].nex = head[u];
    head[u] = e++;
}
  
bool cmp(int i, int j)  {
    return w[i] < w[j];
}
  
void compress(int n)    {
    for (int i=1; i<=n; ++i)    p[i] = i;
    sort (p+1, p+1+n, cmp);
    int rk = 0, pre = w[p[1]] - 1;
    for (int i=1; i<=n; ++i)    {
        if (pre != w[p[i]]) {
            pre = w[p[i]];
            w[p[i]] = ++rk;
        }
        else    {
            w[p[i]] = rk;
        }
    }
}
  
void DFS(int u, int fa) {
    dfn[u] = ++dep; val[dep] = w[u];
    for (int i=head[u]; ~i; i=edge[i].nex)  {
        int v = edge[i].v;
        if (v == fa)    continue;
        DFS (v, u);
    }
    low[u] = dep;
}
  
inline void updata(int x, int c)    {
    cnt[x] += c;
    if (cnt[x] == k)    sum++;
    else if (c > 0 && cnt[x] == k + 1) sum--;
    else if (c < 0 && cnt[x] == k - 1)  sum--;
}
  
void Modui(void)    {
    memset (cnt, 0, sizeof (cnt));
    sum = 0;
    int l = 1, r = 0;
    for (int i=1; i<=m; ++i)    {
        while (data[i].l < l)   updata (val[--l], 1);
        while (data[i].l > l)   updata (val[l], -1), l++;
        while (data[i].r > r)   updata (val[++r], 1);
        while (data[i].r < r)   updata (val[r], -1), r--;
        ans[data[i].id] = sum;
    }
    for (int i=1; i<=m; ++i)    {
        printf ("%d\n", ans[i]);
    }
}
  
int main(void)    {
    int T, cas = 0;  scanf ("%d", &T);
    while (T--) {
        if (cas)    puts ("");
        printf ("Case #%d:\n", ++cas);
        scanf ("%d%d", &n, &k);
        init ();
  
        for (int i=1; i<=n; ++i)    {
            scanf ("%d", &w[i]);
        }
        compress (n);                                           //离散化，升序排序，相同的还是相同的
  
        for (int u, v, i=1; i<n; ++i)   {
            scanf ("%d%d", &u, &v);
            add_edge (u, v);    add_edge (v, u);
        }
        DFS (1, -1);                                            //先序遍历，得到DFS序，树形变线形
  
        int block = (int) sqrt (n + 0.5);
        scanf ("%d", &m);
        for (int u, i=1; i<=m; ++i)    {
            scanf ("%d", &u);
            data[i] = Data (dfn[u], low[u], i, dfn[u] / block);
        }
        sort (data+1, data+1+m);                                //按照DFS序排序
        Modui ();
    }
  
    return 0;
}
/****************************************************************************/

/****************************************************************************/
图论

1 最短路
1.1 Dijkstra
/*
    Dijkstra：单源最短路，优先队列优化，O(ElogV)
    加边用静态链表实现
*/
struct Edge {
    int v, w, nex;
    Edge (int v = 0, int w = 0) : v (v), w (w) {}
    bool operator < (const Edge &r) const {
        return w > r.w;
    }
}edge[E];
bool vis[N];
int head[N], d[N];
int n, m, e;
void Dijkstra(int s)    {
    memset (vis, false, sizeof (vis));
    memset (d, INF, sizeof (d));    d[s] = 0;
    priority_queue<Edge> Q; Q.push (Edge (s, 0));
    while (!Q.empty ()) {
        int u = Q.top ().v;   Q.pop ();
        if (vis[u]) continue;
        vis[u] = true;
        for (int i=head[u]; ~i; i=edge[i].nex)  {
            int v = edge[i].v, w = edge[i].w;
            if (!vis[v] && d[v] > d[u] + w) {
                d[v] = d[u] + w;    Q.push (Edge (v, d[v]));
            }
        }
    }
}
void init(void) {
    memset (head, -1, sizeof (head));
    e = 0;
}
void add_edge(int u, int v, int w)  {
    edge[e].v = v;  edge[e].w = w;
    edge[e].nex = head[u];  head[u] = e++;
}
1.2 SPFA
/*
    SPFA:单源最短路，队列实现，复杂度不定，bellman_Ford 优化版
    return：true：有环；false：无环
*/
bool SPFA(int s)    {
    memset (vis, false, sizeof (vis));
    memset (cnt, 0, sizeof (cnt));
    d[s] = 0;   cnt[s] = 1; vis[s] = true;
    queue<int> Q;   Q.push (s);
    while (!Q.empty ()) {
        int u = Q.front (); Q.pop ();
        vis[u] = false;
        for (int i=head[u]; ~i; i=edge[i].nex)  {
            int v = edge[i].v, w = edge[i].w;
            if (d[v] > d[u] + w) {
                d[v] = d[u] + w;
                if (!vis[v])    {
                    vis[v] = true;  Q.push (v);
                    if (++cnt[v] > n)    return true;
                }
            }
        }
    }
    return false;
}
1.3 Floyd_Warshall
/*
    Floyd_Warshall：任意两点的最短路，动态规划思想，复杂度O(V^3)
    可求传递闭包
*/
void Floyd_Warshall(void)   {
    for (int k=1; k<=n; ++k) {
        for (int i=1; i<=n; ++i) {
            for (int j=1; j<=n; ++j) {
                d[i][j] = min (d[i][j], d[i][k] + d[k][j]);
                // g[i][j] = (g[i][j] || (g[i][k] && g[k][j]));
            }
        }
    }
}

2 拓扑排序
/*
    Topo_Sort：拓扑排序
    in[]：每个点的入度；    ans[]：排序后的结果
    return：true：有环；false：无环
*/
bool Topo_Sort(void)    {
    memset (in, 0, sizeof (in));                                //入度清空
    for (int i=1; i<=n; ++i)    {
        for (int j=0; j<G[i].size (); ++j)  in[G[i][j]]++;      //所有有箭头指向的点的入度累加
    }
    queue<int> Q;   int cnt = 0;
    for (int i=1; i<=n; ++i)    if (!in[i]) Q.push (i);         //先找出入读为0的点，升序入队
    while (!Q.empty ()) {
        int u = Q.front (); Q.pop ();
        ans[++cnt] = u;                                         //ans[]存储拓扑排序后的点
        for (int i=0; i<G[u].size (); ++i)    {
            int v = G[u][i];
            if (!(--in[v])) Q.push (v);                         //每次减少一个入度，若为0，应最先入队，拓扑序靠前
        }
    }
    if (cnt == n)   return false;
    else    return true;
}

3 最小生成树
3.1 Kruskal
/*
    Kruskal：并查集实现，记录两点和距离，按距离升序排序，O (ElogE)
*/
struct Edge {
    int u, v, w;
    bool operator < (const Edge &r) const {
        return w < r.w;
    }
}edge[E];
sort (edge+1, edge+1+m);
if (!uf.same (x, y))    uf.Union (x, y), ans += w;
3.2 Prim
/*
    Prim：Dijkstra思想，优先队列优化，适合稀疏图，O (ElogV)
    不连通返回-1，或返回最小生成树长度(MST)
*/
int Prim(int s) {
    memset (vis, false, sizeof (vis));
    memset (d, INF, sizeof (d));
    priority_queue<Edge> Q;
    for (int i=head[s]; ~i; i=edge[i].nex)  {
        int v = edge[i].v, w = edge[i].w;
        if (d[v] > w)    {
            d[v] = w;   Q.push (Edge (v, d[v]));
        }
    }
    vis[s] = true;  d[s] = 0;   int ret = 0;
    while (!Q.empty ()) {
        int u = Q.top ().v; Q.pop ();
        if (vis[u]) continue;
        vis[u] = true;  ret += d[u];
        for (int i=head[u]; ~i; i=edge[i].nex)  {
            int v = edge[i].v, w = edge[i].w;
            if (!vis[v] && d[v] > w) {
                d[v] = w;   Q.push (Edge (v, d[v]));
            }
        }
    }
    return ret;
}

4 二分图最大匹配
/*
	匈牙利算法,un和vn表示两个集合的点的个数，hungary返回最大匹配数
*/
int un, vn;
int hungary(void)   {
	int res = 0;
    memset (lk, -1, sizeof (lk));
	for (int i=1; i<=un; ++i)   {
		memset (vis, false, sizeof (vis));
		if (DFS (i))	res++;
	}
	return res;
}
bool DFS(int u) {
    for (int i=0; i<G[u].size (); ++i)  {
		int v = G[u][i];
		if (!vis[v])    {
			vis[v] = true;
            if (lk[v] == -1 || DFS (lk[v])) {
				lk[v] = u;	return true;
			}
		}
	}
	return false;
}

5 二分图最大权匹配
/*
	KM算法用来解决最大权匹配问题： 在一个二分图内，左顶点为X，右顶点为Y，现对于每组左右连接Xi，Yj有权w(i,j)，
	求一种匹配使得所有w(i,j)的和最大。也就是最大权匹配一定是完备匹配。如果两边的点数相等则是完美匹配。
	如果点数不相等，其实可以虚拟一些点，使得点数相等，也成为了完美匹配。最大权匹配还可以用最大流去解决
*/
bool DFS(int u)	{		//hungary算法
	visx[u] = true;
	for (int i=1; i<=n; ++i)	{
		if (x[u] + y[i] == w[u][i] && !visy[i])	{
			visy[i] = true;
			if (ly[i] == -1 || DFS (ly[i]))	{
				ly[i] = u;	return true;
			}
		}
		else if (x[u] + y[i] > w[u][i])	d = min (d, x[u] + y[i] - w[u][i]);		//更新d，贪心思想
	}
	return false;
}
int KM(void)	{
	for (int i=1; i<=n; ++i)	{
		x[i] = 0;
		for (int j=1; j<=n; ++j)	{
			x[i] = max (x[i], w[i][j]);			//初始x标杆为最大值w，y为0
		}
	}
	memset (y, 0, sizeof (y));
	memset (ly, -1, sizeof (ly));
	for (int i=1; i<=n; ++i)	{
		while (true)	{
			memset (visx, false, sizeof (visx));
			memset (visy, false, sizeof (visy));
			d = INF;
			if (DFS (i))	break;				//找到增广轨，退出
			for (int j=1; j<=n; ++j)	{		//没有找到，对标杆进行调整
				if (visx[j])	x[j] -= d;
				if (visy[j])	y[j] += d;
			}
		}
	}
	/*
		for (int i=1; i<=n; ++i)    {			//完美匹配的判定
        	if (ly[i] == -1 || w[ly[i]][i] == -INF) return -1;
    	}
    */
	int ret = 0;
	for (int i=1; i<=n; ++i)	{
		ret += x[i] + y[i];
	}
	return ret;
}
/*
	最小费用流：KM算法是求最大流，只要w = -w就可以了，很经典的方法
	即：printf ("%d\n", -KM ());
*/

6 LCA(Lowest Common Ancestor)
/*
    LCA(倍增法，二分搜索)：rt[i][u](i<D=20) 表示u的第2^i的祖先
    LCA预处理复杂度O (logn)，每次询问O (logn)
    DFS中要记录点的深度以及它的父亲，(dep[u] = d; rt[0][u] = fa;)
*/
int LCA(int u, int v)   {
    if (dep[u] < dep[v])    swap (u, v);
    for (int i=0; i<D; ++i) {
        if ((dep[u] - dep[v]) >> i & 1) {
            u = rt[i][u];
        }
    }
    if (u == v) return u;
    for (int i=D-1; i>=0; --i)  {
        if (rt[i][u] != rt[i][v])   {
            u = rt[i][u];
            v = rt[i][v];
        }
    }
    return rt[0][u];
}
void solve(void)    {
    DFS (1, -1, 0);
    for (int i=1; i<D; ++i) {           //init_LCA
        for (int j=1; j<=n; ++j)    {
            rt[i][j] = rt[i-1][j] == -1 ? -1 : rt[i-1][rt[i-1][j]];
        }
    }
    int lca = LCA (u, v);
}
/*
    LCA -> RMQ: LCA (u, v) = F[id[u] <= i <= id[v]中dep[i]最小的i];
    RMQ预处理复杂度O(nlogn)，每次询问O (1)
    (dp[N<<1][D], F[N<<1], dep[N<<1], id[N])
*/
void DFS(int u, int fa, int d, int &k)  {
    id[u] = k;  F[k] = u;   dep[k++] = d;
    for (int i=head[u]; ~i; i=edge[i].nex)  {
        int v = edge[i].v;
        if (v == fa)    continue;
        DFS (v, u, d + 1, k);
        F[k] = u;   dep[k++] = d;
    }
}
int Min(int i, int j)   {
    return dep[i] < dep[j] ? i : j;
}
void init_RMQ(int k)  {
    for (int i=0; i<k; ++i)   dp[i][0] = i;
    for (int j=1; (1<<j)<=k; ++j) {
        for (int i=1; i+(1<<j)-1<k; ++i) {
            dp[i][j] = Min (dp[i][j-1], dp[i+(1<<(j-1))][j-1]);
        }
    }
}
int RMQ(int l, int r)   {
    int k = 0;  while (1<<(k+1) <= r - l + 1)   k++;
    return Min (dp[l][k], dp[r-(1<<k)+1][k]);
}
int LCA(int u, int v)   {
    u = id[u];  v = id[v];
    return u <= v ? F[RMQ (u, v)] : F[RMQ (v, u)];
}
/*
    LCA离线处理，Tarjan算法，复杂度O (N+Q)
    对询问次序按深搜时遍历到的节点顺序进行重组，并查集找祖先
    ans[i]表示第i个询问的LCA
*/
void LCA(int u)    {
    rt[u] = u;
    for (int i=head[u]; ~i; i=edge[i].nex)  {
        int v = edge[i].v;
        if (rt[v] == -1)    {
            LCA (v);
            rt[v] = u;
        }
    }
    for (int i=headq[u]; ~i; i=query[i].nex)    {
        int v = query[i].v;
        if (rt[v] != -1)    {
            int lca = Find (v);
            ans[query[i].id] = lca;
        }
    }
}
/****************************************************************************/

/****************************************************************************/
表达式求值
/*
    表达式求值，逆波兰式(后缀表达式)算法
    输入(可以有空格，支持小数，实现'+-/*%')： ((1+2)*5+1)/4=
    注意：取模一定是要整型，实现版本数字全是double，强制类型转换可能倒置错误
    转换为后缀表达式： 得到：1 2 + 5 * 1 + 4 / =
    计算后缀表达式：得到：4.00
*/
struct Exp  {
    stack<char> op;
    stack<double> num;
    bool error;
 
    int prior(char ch)  {                          //运算符的优先级
        switch (ch) {
            case '+':
            case '-': return 1;
            case '*':
            case '%':
            case '/': return 2;
            default:  return 0;
        }
    }
    bool is_digit(char ch)  {
        return '0' <= ch && ch <= '9';
    }
    string get_postfix(string s)    {              //中缀表达式转变后缀表达式
        while (!op.empty ())    op.pop ();
        op.push ('#');
        string ret = "";
        int len = s.length (), i = 0;
        while (i < len)    {
            if (s[i] == ' ' || s[i] == '=')    {
                i++;    continue;
            }
            else if (s[i] == '(')    {
                op.push (s[i++]);
            }
            else if (s[i] == ')')   {
                while (op.top () != '(')    {
                    ret += op.top ();   ret += ' ';
                    op.pop ();
                }
                op.pop ();  i++;
            }
            else if (s[i] == '+' || s[i] == '-' || s[i] == '*' || s[i] == '/' || s[i] == '%')  {
                while (prior (op.top ()) >= prior (s[i]))    {
                    ret += op.top ();   ret += ' ';
                    op.pop ();
                }
                op.push (s[i++]);
            }
            else    {
                while (is_digit (s[i]) || s[i] == '.')  {
                    ret += s[i++];
                }
                ret += ' ';
            }
        }
        while (op.top () != '#') {
            ret += op.top ();   ret += ' ';
            op.pop ();
        }
        ret += '=';
        return ret;
    }
    double cal(double a, double b, char ch) {
        if (ch == '+')  return a + b;
        if (ch == '-')  return a - b;
        if (ch == '*')  return a * b;
        if (ch == '%')  return (int)((int)a % (int)b);
        if (ch == '/')  {
            if (b != 0) return a / b;
            error = true;   return 0;
        }
    }
    double solve(string str)    {                   //计算后缀表达式
        string s = get_postfix (str);
        while (!num.empty ())   num.pop ();
        error = false;
        int len = s.length (), i = 0;
        while (i < len)  {
            if (s[i] == ' ' || s[i] == '=') {i++;   continue;}
            else if (s[i] == '+' || s[i] == '-' || s[i] == '*' || s[i] == '/' || s[i] == '%')  {
                double a = num.top ();  num.pop ();
                double b = num.top ();  num.pop ();
                num.push (cal (b, a, s[i]));    i++;
            }
            else  {
                double x = 0;
                while (is_digit (s[i])) {
                    x = x * 10 + s[i] - '0';    i++;
                }
                if (s[i] == '.')    {
                    double k = 10.0, y = 0;
                    i++;
                    while (is_digit (s[i])) {
                        y += ((s[i] - '0') / k);
                        i++;    k *= 10;
                    }
                    x += y;
                }
                num.push (x);
            }
        }
        return num.top ();
    }
}E;
int main(void)    {
    ios::sync_with_stdio (false);           //如果全用流的话，加这句话能跑快点
    int T;  cin >> T;
    string str; getline (cin, str);
    while (T--) {
        getline (cin, str);
        cout << E.get_postfix (str) << endl;
        cout << fixed << setprecision (6) << E.solve (str) << endl;
    }
 
    return 0;
}
/****************************************************************************/

/****************************************************************************/
黑科技
1 快速读入输出！
/*
    快速读入输出(读入输出外挂)!--黑科技
    使用场合：huge input (1e6以上)
*/
inline int read(void)   {
    int f = 1, ret = 0; char ch = getchar ();
    while ('0' > ch || ch > '9')    {
        if (ch == '-')  f = -1;
        ch = getchar ();
    }
    while ('0' <= ch && ch <= '9')  {
        ret = ret * 10 + ch - '0';
        ch = getchar ();
    }
    return ret * f;
}
inline void print(int x)    {
    if (x < 0)  {
        putchar ('-');  x = -x;
    }
    if (x > 9)  print (x / 10);
    putchar (x % 10 + '0');
}

2 解决爆栈，手动加栈！
#pragma comment(linker, "/STACK:1024000000,1024000000")

3 #include <bits/stdc++.h>
/****************************************************************************/
