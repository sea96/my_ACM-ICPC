/*
    快速读入输出(读入输出外挂)!--黑科技
    使用场合：huge input (1e6以上)
*/
inline int read() {
    int f = 1, ret = 0;
    char ch = getchar ();
    while ('0' > ch || ch > '9') {
        if (ch == '-')  f = -1;
        ch = getchar ();
    }
    while ('0' <= ch && ch <= '9') {
        ret = ret * 10 + ch - '0';
        ch = getchar ();
    }
    return ret * f;
}
inline void print(int x) {
    if (x < 0) {
        putchar ('-');
        x = -x;
    }
    if (x > 9) {
        print (x / 10);
    }
    putchar (x % 10 + '0');
}
