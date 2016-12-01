//辗转相除法(欧几里德算法)
int GCD(int a, int b) {
    return b ? GCD (b, a % b) : a;
}
//快速GCD
int quick_GCD(int a, int b) {
    if (!a) return b;
    if (!b) return a;
    if (! (a & 1) && ! (b & 1)) return quick_GCD (a >> 1, b >> 1) << 1;
    else if (! (b & 1)) return quick_GCD (a, b >> 1);
    else if (! (a & 1)) return quick_GCD (a >> 1, b);
    else    return quick_GCD (abs (a - b), min (a, b));
}

/*
 *扩展欧几里得算法，求x, y，使得 ax + by = GCD(a, b)
 *返回d = GCD(a,b)，和对应于等式 ax + by = d中的x，y
 *求解 ax + by = c: 当c是d的倍数时，记倍数为k，解为x0=kx，y0=ky；否则无解
 *通解：x = x0 + (b/d)*n, y = y0 - (a/d)*n, n为整数
 */
void ex_GCD(int a, int b, int &x, int &y, int &d) {
    if (!b) { x = 1; y = 0; d = a; }
    else {
        ex_GCD (b, a % b, y, x, d);
        y -= a / b * x;
    }
}

/*
 *模线性方程(modular_linear_equation) ax = b (mod n)
 *转换为ax = ny + b，即ax + ny = b。答案保存在res里
 *返回true，有解；false，无解
 */
bool linear_mod_equ(int a, int b, int mod, vector<int> &res) {
    int x, y, d;
    ex_GCD(a, mod, x, y, d);
    if (b % d) return false;
    x = x * (b/d) % mod;
    res.push_back(x);
    for (int i=1; i<d; ++i) {
        res.push_back((x + (int)i*(mod/d)) % mod);
    }
    return true;
}
