//快速读入，输出
inline void read(int &x) {
    char ch = getchar ();
    int f = 1; x = 0;
    while (ch < '0' || ch > '9') if (ch == '-') f = -1, ch = getchar ();
    while (ch >= '0' && ch <= '9') x = x * 10 + (ch - '0'), ch = getchar ();
    x *= f;
}
inline void write(int x) {
    if (x < 0) putchar ('-'), x = -x;
    if (x > 9) write (x / 10);
    putchar (x % 10 + '0');
}
