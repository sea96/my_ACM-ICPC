//递推，对MOD没有要求，预处理时间复杂度O(n^2)
//C[n][m] = C (n, m) % MOD
void init() {
    for (int i=0; i<N; ++i) {
        C[i][0] = C[i][i] = 1;
        for (int j=1; j<i; ++j) {
            C[i][j] = C[i-1][j] + C[i-1][j-1];
            if (C[i][j] >= MOD) {
                C[i][j] -= MOD;
            }
        }
    }
}

//C(n, k)，分子和分母分别累乘，时间复杂度O(min(k, n-k))
long long binom(long long n, long long k) {
    if (k > n - k) k = n - k;
    long long a = 1, b = 1;
    for (int i=1; i<=k; ++i) {
        a *= (n + 1 - i);
        b *= i;
        if (a % b == 0) { a /= b; b = 1; }
    }
    return a / b;
}

//要求MOD是质数，预处理时间复杂度O(n)
const int N = 1e5 + 10;
const int MOD = 1e9 + 7;
int f[N], finv[N], inv[N];
void init() {
    inv[1] = 1;
    for (int i=2; i<N; ++i) {
        inv[i] = (MOD - MOD / i) * 1ll * inv[MOD%i] % MOD;
    }
    f[0] = finv[0] = 1;
    for (int i=1; i<N; ++i) {
        f[i] = f[i-1] * 1ll * i % MOD;
        finv[i] = finv[i-1] * 1ll * inv[i] % MOD;
    }
}

//C(n, k) % p
ll comb(ll n, ll k, ll p)  {
    if (k < 0 || k > n) return 0;
    return f[n] * inv[f[n-k]*f[k]%p] % p;
}
//Lucas定理，大组合数取模。要求p是质数
//如果p不是质数，质因数分解后用Lucas定理，然后用中国剩余定理（CRT）。
ll Lucas(ll n, ll k, ll p) {
    if (n < k) return 0;
    if (!k) return 1;
    return comb(n%p, k%p, p, id) * Lucas(n/p, k/p, p) % p;
}
