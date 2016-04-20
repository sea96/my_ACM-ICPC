//辗转相除法(欧几里德算法)
int GCD(int a, int b) {
    return b ? GCD (b, a % b) : a;
}
//非递归版本
int GCD(int a, int b) {
    while (b) {
        int c = b;
        b = a % b;
        a = c;
    }
    return a;
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
    *扩展欧几里得算法，求x, y 使得 ax + by = GCD (a, b)
    *返回d = GCD (a,b)，和对应于等式 ax + by = d中的x，y
    *求解 ax + by = c: 当c是d的倍数时，记倍数为k，解为kx，ky；否则无解
*/
int ex_GCD(int a, int b, int &x, int &y)    {
    if (!a && !b)   return -1;  //无最大公约数
    if (!b) {
        x = 1; y = 0;   //a*1+0*b=a
        return a;
    }
    int d = ex_GCD (b, a % b, y, x);
    y -= a / b * x;
    return d;
}
/*
    *模线性方程 ax = b (mod n)
    *转换为ax = ny + b，即ax - ny = b。答案保存在xs里
    *返回true，有解；false，无解
*/
bool modular_linear_equation(int a, int b, int n, std::vector<int> &xs) {
    int x, y;
    int d = ex_GCD (a, n, x, y);
    if (b % d) {
        return false;
    }
    x = x * (b / d) % n;
    for (int i=1; i<=d; ++i) {
        xs.push_back ((x + i * (n/d)) %n);
    }
    return true;
}
