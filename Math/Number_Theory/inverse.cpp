//求逆元，三种方法都要求mod为质数
//求法1：预处理，数组保存，复杂度 O(n)
void Inv(int n, int mod) {
    inv[1] = 1;
    for (int i=2; i<=n; ++i) {
        inv[i] = 1ll * (mod - mod / i) * inv[mod%i] % mod;
    }
}
//求法2：复杂度 O(logn)
int Inv(int a, int mod) {
    return pow_mod (a, mod - 2, mod);
}
//求法3：扩展欧几里德求解
int Inv(int a, int mod) {
    int x, y;
    int d = ex_GCD (a, mod, x, y);
    if (d == 1) {
        return (x % mod + mod) % mod;
    } else {
        return -1;
    }
}
