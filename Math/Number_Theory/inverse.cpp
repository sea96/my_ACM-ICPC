//求逆元，三种方法都要求p为质数
//求法1：预处理，数组保存，复杂度 O (n)
void Inv(int n, int p) {
    inv[1] = 1;
    for (int i=2; i<=n; ++i) {
        inv[i] = 1ll * (p - p / i) * inv[p%i] % p;
    }
}
//求法2：复杂度 O (logn)
int Inv(int a, int p) {
    return pow_mod (a, p - 2, p);
}
//求法3：扩展欧几里德求解
int Inv(int a, int p) {
    int x, y;
    int d = ex_GCD (a, p, x, y);
    if (d == 1) return (x % p + p) % p;
    else    return -1;
}
