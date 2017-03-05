//快速读入(int,long long,double)
template<class T>
inline void read(T &x) {
    x = 0; static char c; bool minus = false;
    for (; !(c>='0' && c<='9'); c=getchar()) if (c == '-') minus = true;
    for (; c>='0' && c<='9'; x = x*10+(c-'0'), c=getchar());
    if (c == '.') {
        double tmp = 1.0;
        c = getchar();
        for (; c>='0' && c<='9'; tmp /= 10, x += tmp*(c-'0'), c=getchar());
    }
    if (minus) x = -x;
}

//快速读入（fread版），有对EOF的处理
inline char getc() {
    static char buf[1000000], *p1 = buf, *p2 = buf;
    if (p1 == p2) {
        p2 = (p1 = buf) + fread(buf, 1, 1000000, stdin);
        if (p1 == p2) return EOF;
    }
    return *p1++;
}
template<class T>
inline void read(T &x) {
    x = 0; static char c; bool minus = false;
    for (; !(c>='0' && c<='9'); c=getc()) {
        if (c == '-') minus = true;
        if (c == EOF) { x = INT_MIN; return ; }
    }
    for (; c>='0' && c<='9'; x = x*10+(c-'0'), c=getc()); if (minus) x = -x;
}

//快速读入(char*)
inline bool blank(char ch) {
    return ch==' ' || ch=='\n' || ch=='\r' || ch=='\t';
}
inline void read(char *s) {
    char c = getchar();
    for (; blank(c); c=getchar());
    for (; !blank(c); c=getchar()) *s++ = c;
    *s = 0;
}

//快速输出
inline void write(int x) {
    if (x < 0) putchar ('-'), x = -x;
    if (x > 9) write (x / 10);
    putchar (x % 10 + '0');
}
