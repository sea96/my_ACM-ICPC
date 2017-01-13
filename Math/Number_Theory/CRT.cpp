/*
 *中国剩余定理（Chinese_remainder_theorem）：最小的x，x=a[i](mod m[i]) (1<=i<=n)
 *做法1（互质）：令M=m[i]的乘积，w[i]=M/m[i]，GCD(w[i], m[i])=1，求解Inv(w[i])（w[i]*x+m[i]*y=1）
 *则ans = sigma(a[i]*w[i]*x);
 *做法2（不互质）：每次联立两个方程组合并得到新的方程组，重复n-1次
 */
ll CRT(ll *a, ll *m, int n) {
    ll M = 1, x, y, z, ret = 0;
    for (int i=1; i<=n; ++i) M *= m[i];
    for (int i=1; i<=n; ++i) {
        ll w = M / m[i];
        ex_GCD(w, m[i], x, y, z);
        ret = (ret + a[i]*w%M*x%M) % M;
    }
    return (ret + M) % M;
}
ll ex_CRT(ll *a, ll *m, int n) {
    ll lcm = m[0], r = a[0], x, y, d;
    for (int i=1; i<n; ++i) {
        ex_GCD(lcm, m[i], x, y, d);
        if ((a[i]-r) % d) return -1;
        x = (a[i]-r) / d * x % (m[i]/d);
        r += x * lcm;
        lcm = lcm / d * m[i];
        r %= lcm;
    }
    return (r + lcm) % lcm;
}
