//欧拉筛法(Euler)：返回n以内素数的个数（保存在prime[0]）， 复杂度O(n) ?
void prime_table(int n) {
    memset (is_prime, true, sizeof (is_prime));
    is_prime[0] = is_prime[1] = false;
    int &p = prime[0] = 0;
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
 *素性测试，Miller_Rabin 测试算法
 *if n < 3,215,031,751, just test a = 2, 3, 5, and 7;
 *if n < 4,759,123,141, just test a = 2, 7, and 61;
 *if n < 1,122,004,669,633, just test a = 2, 13, 23, and 1662803;
 *if n < 2,152,302,898,747, just test a = 2, 3, 5, 7, and 11;
 *if n < 3,474,749,660,383, just test a = 2, 3, 5, 7, 11, and 13;
 *if n < 341,550,071,728,321, just test a = 2, 3, 5, 7, 11, 13, and 17.
 *if n < 3,825,123,056,546,413,051, just test a = 2, 3, 5, 7, 11, 13, 17, 19, and 23.
 *if n < 18,446,744,073,709,551,616 = 2^64, just test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
 */
bool test(ll x, ll n) {
    ll m = n - 1;
    while (!(m&1)) m >>= 1;
    x = pow_mod(x, m, n);
    for (; m<n-1 && x!=1 && x!=n-1; m<<=1) {
        x = mul_mod(x, x, n);
    }
    return x == n-1 || (m&1) == 1;
}
bool Miller_Rabin(ll n) {
    if (n == 2 || n == 7 || n == 61) return true;
    if (n < 2 || !(n&1)) return false;
    return test(2, n) && test(7, n) && test(61, n);
}

//唯一分解定理，先打个素数表优化试除法
void get_factors(int n) {
    std::vector<Node> factors;
    for (int i=1; i<=prime[0]; ++i) {
        if (n % prime[i] == 0) {
            int f = prime[i], c = 0;
            while (n % prime[i] == 0) {
                n /= prime[i];
                c++;
            }
            factors.push_back ((Node) {f, c});
        }
    }
    if (n > 1) {
        factors.push_back ((Node) {(n, 1});
    }
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
void factorize(ll n, vector<ll> &ret) {  //ret保存质因数，无序
    if (Miller_Rabin (n)) {  //素数
        ret.push_back(n);  return ;
    }
    ll p = n;
    while (p >= n)   p = Pollard_rho(p, rand () % (n - 1) + 1);
    factorize(p, ret);
    factorize(n / p, ret);
}

//约数枚举，复杂度O (sqrt(n)) By TiaoZhan
vector<int> divisor(int n)  {
    vector<int> res;
    for (int i=1; i*i<=n; ++i)  {
        if (n % i == 0) {
            res.push_back(i);
            if (n / i != i) res.push_back(n / i);
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
//单独求数 欧拉函数phi[n]
int phi(int n) {
    int m = (int)sqrt (n + 0.5);
    int ret = n;
    for (int i=2; i<=m; ++i) {
        if (n % i == 0) {
            ret = ret / i * (i - 1);
            while (n % i == 0) {
                n /= i;
            }
        }
    }
    if (n > 1) {
        ret = ret / n * (n - 1);
    }
    return ret;
}
//线性筛，复杂度 O(nloglogn)，还有欧拉筛素数顺便求phi[i]
void phi_table(int n) {
    memset(phi, 0, sizeof (phi));
    phi[1] = 1;
    for (int i=2; i<=n; ++i) {
        if (!phi[i]) {
            for (int j=i; j<=n; j+=i) {
                if (!phi[j]) {
                    phi[j] = j;
                }
                phi[j] = phi[j] / i * (i - 1);
            }
        }
    }
}
