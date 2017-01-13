/*
 *中国剩余定理（Chinese_remainder_theorem）：最小的x，x=a[i](mod m[i]) (1<=i<=k)
 */
ll CRT(ll *a, ll *m, int k) {
    ll M = 1, x, y, z, ret = 0;
    for (int i=1; i<=k; ++i) M *= m[i];
    for (int i=1; i<=k; ++i) {
        ll w = M / m[i];
        ex_GCD(w, m[i], x, y, z);
        ret = (ret + x*w*a[i]) % M;
    }
    return (ret + M) % M;
}
/*
 *做法：令M=m[i]的乘积，w[i]=M/m[i]，GCD(w[i], m[i])=1，求解w[i]*x+m[i]*y=1
 *令e[i]=w[i]*x，ans = sum(e[i]*a[i]);
 *求解满足：x = 2(% 3), x = 3(% 5), x = 2(% 7)的最小x
 *a对应(2，3，2)，m对应(3，5，7)
 *解为：x = 2*35a + 3*21b + 2*15c (%105)
 */
