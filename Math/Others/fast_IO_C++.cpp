//快速读入，输出
inline void read(int &x) {
    char ch = getchar ();
    bool sign = false; x = 0;
    for (; ch < '0' || ch > '9'; ch=getchar ()) if (ch == '-') sign = true;
    for (; ch >= '0' && ch <= '9'; ch=getchar ()) x = x * 10 + (ch - '0');
    if (sign) x = -x;
}

inline void read(double &x) {
    char ch = getchar ();
    bool sign = false; x = 0;
    for (; ch < '0' || ch > '9'; ch=getchar ()) if (ch == '-') sign = true;
    for (; ch >= '0' && ch <= '9'; ch=getchar ()) x = x * 10 + (ch - '0');
    if (ch == '.') {
        double tmp = 1;
        ch = getchar ();
        for (; ch >= '0' && ch <= '9'; ch=getchar ()) tmp /= 10, x += tmp * (ch - '0');
    }
    if (sign) x = -x;
}

inline bool blank(char ch) {
    return ch==' ' || ch=='\n' || ch=='\r' || ch=='\t';
}

inline void read(char *s) {
    char ch = getchar ();
    for (; blank (ch); ch=getchar ());
    for (; !blank (ch); ch=getchar ()) *s++ = ch;
    *s = 0;
}

inline void write(int x) {
    if (x < 0) putchar ('-'), x = -x;
    if (x > 9) write (x / 10);
    putchar (x % 10 + '0');
}
