//递推，对MOD没有要求，预处理时间复杂度O(n^2)
//C[n][m] = C (n, m) % MOD
void init() {
    C[0][0] = 1;
    for (int i=0; i<25; ++i) {
        C[i][0] = C[i][i] = 1;
        for (int j=1; j<i; ++i) {
            C[i][j] = C[i-1][j] + C[i-1][j-1];
            if (C[i][j] >= MOD) {
                C[i][j] -= MOD;
            }
        }
    }
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
//C (n, k) % MOD
int comb(int n, int k)  {
    if (k < 0 || k > n) return 0;
    return f[n] * 1ll * finv[n-k] % MOD * finv[k] % MOD;
}

//Lucas定理，大组合数取模
ll f[N];　　//f[n] = n!
void init(int p)　{
    f[0] = 1;
    for (int i=1; i<=p; ++i) f[i] = f[i-1] * i % p;
}
//C (n, k) % p
ll Lucas(ll n, ll k, int p)　{
     ll ret = 1;
     while (n && k) {
        ll nn = n % p, kk = k % p;
        if (nn < kk) return 0;                   //inv (f[kk]) = f[kk] ^ (p - 2) % p
        ret = ret * f[nn] * pow_mod (f[kk] * f[nn-kk] % p, p - 2, p) % p;
        n /= p, k /= p;
     }
     return ret;
}
