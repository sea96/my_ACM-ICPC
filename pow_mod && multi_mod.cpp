//乘方取模(快速幂)，x ^ n % p，复杂度O (logn)
int pow_mod(int x, int n, int p) {
    int ret = 1;
    while (n) {
        if (n & 1) {
            ret = 1ll * ret * x % p;
        }
        x = 1ll * x * x % p;
        n >>= 1;
    }
    return ret;
}
//乘法取模，a * b % p，复杂度O (logb)
int multi_mod(int a, int b, int p) {
    int ret = 0;
    a = (a % p + p) % p;
    b = (b % p + p) % p;
    while (b) {
        if (b & 1) {
            ret += a;
            if (ret >= p)   ret -= p;
        }
        b >>= 1;
        a <<= 1;
        if (a >= p) a -= p;
    }
    return ret;
}
