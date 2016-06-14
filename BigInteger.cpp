/*
    *大数模版 copy from JayYe
*/
const int numlen = 2005; // 需要的位数
const int numbit = 4;   // 数组一位表示的整数
const int addbit = 10000;//进位数
const int maxn = numlen/numbit + 10;   // 数组要开的位数

int Max(int a, int b) { return a>b?a:b; }
struct bign {
    int len, s[numlen];
    bign() {
        memset(s, 0, sizeof(s));
        len = 1;
    }
    bign(int num) { *this = num; }
    bign(const char *num) { *this = num; }
    bign operator = (const int num) {
        char s[numlen];
        sprintf(s, "%d", num);
        *this = s;
        return *this;
    }
    bign operator = (const char *num) {
        int clen = strlen(num);
        while(clen > 1 && num[0] == '0') num++, clen--;
        len = 0;
        for(int i = clen-1;i >= 0;i -= numbit) {
            int top = min(numbit, i+1), mul = 1;
            s[len] = 0;
            for(int j = 0;j < top; j++) {
                s[len] += (num[i-j]-'0')*mul;
                mul *= 10;
            }
            len++;
        }
        deal();
        return *this;
    }

    void deal() {
        while(len > 1 && !s[len-1]) len--;
    }

    bign operator + (const bign &a) const {
        bign ret;
        ret.len = 0;
        int top = Max(len, a.len) , add = 0;
        for(int i = 0;add || i < top; i++) {
            int now = add;
            if(i < len) now += s[i];
            if(i < a.len)   now += a.s[i];
            ret.s[ret.len++] = now%addbit;
            add = now/addbit;
        }
        return ret;
    }
    bign operator - (const bign &a) const {
        bign ret;
        ret.len = 0;
        int cal = 0;
        for(int i = 0;i < len; i++) {
            int now = s[i] - cal;
            if(i < a.len)   now -= a.s[i];
            if(now >= 0)    cal = 0;
            else {
                cal = 1; now += addbit;
            }
            ret.s[ret.len++] = now;
        }
        ret.deal();
        return ret;
    }
    bign operator * (const bign &a) const {
        bign ret;
        ret.len = len + a.len;
        for(int i = 0;i < len; i++) {
            int pre = 0;
            for(int j = 0;j < a.len; j++) {
                int now = s[i]*a.s[j] + pre;
                pre = 0;
                ret.s[i+j] += now;
                if(ret.s[i+j] >= addbit) {
                    pre = ret.s[i+j]/addbit;
                    ret.s[i+j] -= pre*addbit;
                }
            }
            if(pre) ret.s[i+a.len] = pre;
        }
        ret.deal();
        return ret;
    }
    // 除以一个小整数     ***********注意计算过程可能会爆int
    bign operator / (const int a) const {
        bign ret;
        ret.len = len;
        int pre = 0;
        for(int i = len-1;i >= 0; i--) {
            ret.s[i] = (s[i] + pre*addbit)/a;
            pre = s[i] + pre*addbit - a*ret.s[i];
        }
        ret.deal();
        return ret;
    }
    bign operator % (const int a) const {
        bign b = *this / a;
        return *this - b*a;
    }

    bign operator += (const bign &a) { *this = *this + a; return *this; }
    bign operator -= (const bign &a) { *this = *this - a; return *this; }
    bign operator *= (const bign &a) { *this = *this * a; return *this; }
    bign operator /= (const int a) { *this = *this / a; return *this; }
    bign operator %= (const int a) { *this = *this % a; return *this; }

    bool operator < (const bign &a) const {
        if(len != a.len)    return len < a.len;
        for(int i = len-1;i >= 0; i--) if(s[i] != a.s[i])
            return s[i] < a.s[i];
        return false;
    }
    bool operator > (const bign &a) const  { return a < *this; }
    bool operator <= (const bign &a) const { return !(*this > a); }
    bool operator >= (const bign &a) const { return !(*this < a); }
    bool operator == (const bign &a) const { return !(*this > a || *this < a); }
    bool operator != (const bign &a) const { return *this > a || *this < a; }
};
istream& operator >> (istream &in, bign &x) {
    string s;
    in >> s;
    x = s.c_str();
    return in;
}
ostream& operator << (ostream &out, const bign &x) {
    printf("%d", x.s[x.len-1]);
    for(int i = x.len-2;i >= 0; i--)    printf("%04d", x.s[i]);
    return out;
}

bign f[110];    //卡特兰数
int main(void)  {
    f[0] = 1;
    for (int i=1; i<=100; ++i)  {
        f[i] = f[i-1] * (4 * i - 2) / (i + 1);
    }
    int n;
    while (scanf ("%d", &n) == 1)   {
        if (n == -1)    break;
        cout << f[n] << endl;
    }

    return 0;
}

/*
    大数介乘求位数
    log10(n!) = log10(1*2*..*n) = log10(1) + log10(2) + ...+log10(n)
    解释：123456=1.23456*10^5;
          log10(123456)=5.09151;
          log10(1.23456*10^5)=log10(1.23456)+log10(10^5)=0.09151+5;
          故int(log10(i))+1 (i<=n) 就是n!的位数
*/
