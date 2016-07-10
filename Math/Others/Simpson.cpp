double F(double x) {
    //待积分函数
}
//三点Simpson法，要求F是全局函数
double Simpson(double a, double b) {
    double c = a + (b - a) / 2;
    return (F (a)+4*F (c)+F (b)) * (b-a) / 6;
}
//自适应Simpson公式（递归过程）。整个区间[a,b]的三点Simpson值A
double asr(double a, double b, double eps, double A) {
    double c = a + (b - a) / 2;
    double L = Simpson (a, c), R = Simpson (c, b);
    if (std::fabs (L+R-A) <= 15*eps) {
        return L+R+(L+R-A)/15.0;
    }
    return asr (a, c, eps/2, L) + asr (c, b, eps/2, R);
}
//自适应Simpson公式（主过程）
double asr(double a, double b, double eps) {
    return asr (a, b, eps, Simpson (a, b));
}
