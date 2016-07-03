//欧拉筛法(Euler)：返回n以内素数的个数（保存在prime[0]）， 复杂度O(n) ?
void sieve(int n) {
    int &p = prime[0] = 0;
    memset (is_prime, true, sizeof (is_prime));
    is_prime[0] = is_prime[1] = false;
    //phi[1] = 1;
    for (int i=2; i<=n; ++i) {
        if (is_prime[i]) {
            prime[++p] = i;
            //phi[i] = i - 1;
        }
        for (int j=1; j<=p && i*prime[j]<=n; ++j) {
            is_prime[i*prime[j]] = false;
            if (i % prime[j] == 0) {
                //phi[i*prime[j]] = phi[i] * prime[j];
                break;
            } else {
                //phi[i*prime[j]] = phi[i] * (prime[j] - 1);
            }
        }
    }
}
/*
    *素性测试，在小范围(1e5)内判素数个数以及单个数判素数有奇效
    *不适用于大范围判素数
*/
bool is_prime(int x) {
    if (x == 2 || x == 3)   return true;
    if (x % 6 != 1 && x % 6 != 5)   return false;
    for (int i=5; i*i<=x; i+=6) {
        if (x % i == 0 || x % (i + 2) == 0) return false;
    }
    return x != 1;
}
/*
    *素性测试，对sqrt(n)内的素数进行测试即可，预处理求出sqrt(n)中的素数，
    *假设该范围内素数的个数为s，那么复杂度降为O(s)。
*/
bool is_prime(int x)   {
    int p = seive (sqrt (x + 0.5));
    for (int i=1; i<=p; ++i) {
        if (x % prime[i] == 0) {
            return false;
        }
    }
    return x != 1;
}
/*
    *素性测试，Miller_Rabin 随机算法
    *可以判断< 2^63的数
    以a为基,n-1=x*2^t，a^(n-1) = 1(mod n)  验证n是不是素数
    *返回true：素数；false：合数
*/
bool Miller_Rabin(ll n) {
    if (n == 2) return true;
    if (n < 2 || ! (n & 1))  return false;           //偶数或1
    ll x = n - 1;   int t = 0;
    while (! (x & 1)) {
        x >>= 1;  t++;
    }
    for (int i=1; i<=S; ++i) {                      //srand (time (NULL)); S=20
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

//整数分解，先打个素数表优化试除法
vector<int> factorize(int n) {
    int p = seive (n + 0.5);
    vector<int> ret;
    for (int i=1; i<=p; ++i) {
        while (n % prime[i] == 0) {
            n /= prime[i];
            ret.push_back (prime[i]);
        }
    }
    if (n != 1) ret.push_back (n);
    return ret;
}
/*
    *大整数分解，Pollard_rho 随机算法
    *factorize ()保存质因数在vector
*/
ll Pollard_rho(ll x, ll c) {
    ll i = 1, k = 2;
    ll a = rand () % x;
    ll b = a;
    while (1) {
        i++;
        a = (multi_mod (a, a, x) + c) % x;
        ll d = GCD (b - a, x);
        if (d != 1 && d != x)   return d;
        if (b == a) return x;
        if (i == k) b = a, k += k;
    }
}
void factorize(ll n, vector<ll> &ret) {           //ret保存质因数，无序
    if (Miller_Rabin (n)) {                       //素数
        ret.push_back (n);  return ;
    }
    ll p = n;
    while (p >= n)   p = Pollard_rho (p, rand () % (n - 1) + 1);
    factorize (p, ret);
    factorize (n / p, ret);
}

//约数枚举，复杂度O (sqrt(n)) By TiaoZhan
vector<int> divisor(int n)  {
    vector<int> res;
    for (int i=1; i*i<=n; ++i)  {
        if (n % i == 0) {
            res.push_back (i);
            if (n / i != i) res.push_back (n / i);
        }
    }
    return res;
}

/*
    *对正整数n，欧拉函数是少于或等于n的数中与n互质的数的数目。
    *例如euler(8)=4，因为1,3,5,7均和8互质。
    *Euler函数表达通式：euler(x)=x(1-1/p1)(1-1/p2)(1-1/p3)(1-1/p4)…(1-1/pn),
    *其中p1,p2……pn为x的所有素因数，x是不为0的整数。euler(1)=1 
    *欧拉公式的延伸：一个数的所有质因子之和是euler(n)*n/2。
*/
//线性筛 详见欧拉筛素数，欧拉函数 phi[i]
//单独求数 欧拉函数 phi[n]
int eular(int n)  {
    int ret = n;
    for (int i = 2; i*i<=n; ++i) {
        if (n % i == 0) {
            ret -= ret / i;
            while (n % i == 0) {
                n /= i;
            }
        }
    }
    if (x > 1) ret -= ret / n;
    return ret;
}
