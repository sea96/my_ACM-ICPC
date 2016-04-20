/*
    *中国剩余定理：求解满足：x = 2 (%3), x = 3 (%5), x = 2 (%7)的最小x
    *解为：x = 2 * 35a + 3 * 21b + 2 * 15c (%105)
    *k是方程组的个数，b对应(2，3，2)，m对应(3，5，7)
*/
int China(int k, int *b, int *m) {
    ll M = 1, x, y, ret = 0;
    for (int i=1; i<=k; ++i) {
        M *= m[i];
    }
    for (int i=1; i<=k; ++i) {
        ll w = M / m[i];
        ex_GCD (w, m[i], x, y);
        ret += multi_mod (multi_mod (x, w, M), b[i], M); //防爆long long
    }
    return (ret % M + M) % M;
}
