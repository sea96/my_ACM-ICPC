/*
 *支持大数加减乘除（取余）四则运算，支持小数加减乘除（取余）运算
 *支持大数之间比较大小，（大数与小数的比较将小数转大数）
 */
struct BigInteger {
    static const int BASE = 100000000;  //一亿进制（10^8），即每位0~99999999的数
    static const int DIGIT = 8;  //每位的长度
    static const int LEN = 1000;  //位数，即最多能保存LEN*DIGIT的长度

    int s[LEN];
    int len;

    //构造函数
    BigInteger() { memset (s, 0, sizeof (s)); len = 1; }
    BigInteger(const int rhs) { (*this) = rhs; }
    BigInteger(const char *str) { (*this) = str; }
    BigInteger operator = (int rhs) {
        for (len=0; rhs; rhs/=BASE) s[++len] = rhs % BASE;
        if (!len) s[++len] = 0;
        return *this;
    }
    BigInteger operator = (const char *str) {
        int j = strlen (str) - 1;
        len = j / DIGIT + 1;
        for (int i=0; i<=len; ++i) s[i] = 0;
        for (int i=0; i<=j; ++i) {
            int k = (j - i) / DIGIT + 1;
            s[k] = s[k] * 10 + str[i] - '0';
        }
        return *this;
    }

    //大数与大数的四则运算
    BigInteger operator + (const BigInteger &rhs) const {
        BigInteger ret;
        int i;
        for (i=1; i<=len || i<=rhs.len || ret.s[i]; ++i) {
            if (i <= len) ret.s[i] += s[i];
            if (i <= rhs.len) ret.s[i] += rhs.s[i];
            ret.s[i+1] = ret.s[i] / BASE;
            ret.s[i] %= BASE;
        }
        ret.len = i - 1;
        if (ret.len == 0) ret.len = 1;
        return ret;
    }
    BigInteger operator - (const BigInteger &rhs) const {
        BigInteger ret;
        for (int i=1, j=0; i<=len; ++i) {
            ret.s[i] = s[i] - j;
            if (i <= rhs.len) ret.s[i] -= rhs.s[i];
            if (ret.s[i] < 0) { j = 1; ret.s[i] += BASE; }
            else j = 0;
        }
        ret.len = len;
        while (ret.len > 1 && !ret.s[ret.len]) ret.len--;
        return ret;
    }
    BigInteger operator * (const BigInteger &rhs) const {
        BigInteger ret;
        ret.len = len + rhs.len;
        ll g = 0;
        for (int k=1; k<=ret.len; ++k) {
            ll tmp = g;
            int i = k + 1 - rhs.len;
            if (i < 1) i = 1;
            for (; i<=k && i<=len; ++i)
                tmp += (ll) s[i] * rhs.s[k+1-i];
            g = (int) (tmp / BASE);
            ret.s[k] = (int) (tmp % BASE);
        }
        while (ret.len > 1 && !ret.s[ret.len]) ret.len--;
        return ret;
    }
    //大数除法，返回<商,余数>
    pair<BigInteger, BigInteger> operator / (const BigInteger &rhs) const {
        BigInteger a, b;
        for (int i=len; i>=1; --i) {
            a.s[i] = 0;
            b = b * BASE; b.s[1] = s[i];
            int l = 0, r = BASE - 1, mid;
            while (l < r) {
                mid = (l + r + 1) >> 1;
                if (rhs * mid <= b) l = mid;
                else r = mid - 1;
            }
            a.s[i] = l; b = b - rhs*l;
        }
        a.len = len;
        while (a.len > 1 && !a.s[a.len]) a.len--;
        return make_pair (a, b);
    }

    //大数与小数的四则运算
    BigInteger operator + (const int rhs) const {BigInteger tmp; tmp=rhs; return *this+tmp; }
    BigInteger operator - (const int rhs) const {BigInteger tmp; tmp=rhs; return *this-tmp; }
    BigInteger operator * (const int rhs) const {BigInteger tmp; tmp=rhs; return *this*tmp; }
    BigInteger operator / (const int rhs) const {
        ll d = 0;
        BigInteger ret;
        for (int i=len; i>=1; --i) {
            d = d * BASE + s[i];
            ret.s[i] = d / rhs; d %= rhs;
        }
        ret.len = len;
        while (ret.len > 1 && !ret.s[ret.len]) ret.len--;
        return ret;
    }
    int operator % (const int rhs) const {
        ll d = 0;
        for (int i=len; i>=1; --i) d = (d*BASE+s[i]) % rhs;
        return (int) d;
    }

    BigInteger operator += (const BigInteger &rhs) { *this = *this+rhs; return *this; };
    BigInteger operator -= (const BigInteger &rhs) { *this = *this-rhs; return *this; };
    BigInteger operator *= (const BigInteger &rhs) { *this = *this*rhs; return *this; };
    BigInteger operator /= (const BigInteger &rhs) { *this = (*this/rhs).first; return *this; };
    BigInteger operator %= (const BigInteger &rhs) { *this = (*this/rhs).second; return *this; };

    bool operator < (const BigInteger &rhs) const {
        if (len != rhs.len) return len < rhs.len;
        for (int i=len; i>=1; --i) {
            if (s[i] != rhs.s[i]) return s[i] < rhs.s[i];
        }
        return false;
    }
    bool operator > (const BigInteger &rhs) const { return rhs < *this; }
    bool operator == (const BigInteger &rhs) const { return !(*this < rhs || *this > rhs); }
    bool operator <= (const BigInteger &rhs) const { return !(*this > rhs); }
    bool operator >= (const BigInteger &rhs) const { return !(*this < rhs); }
    bool operator != (const BigInteger &rhs) const { return !(*this == rhs); }
};

istream& operator >> (istream &in, BigInteger &x) {
    string str;
    if (!(in >> str)) return in;
    x = str.c_str ();
    return in;
}

ostream& operator << (ostream &out, const BigInteger &x) {
    printf ("%d", x.s[x.len]);
    for (int i=x.len-1; i>=1; --i) printf ("%08d", x.s[i]);
    return out;
}

BigInteger f[105];  //卡特兰数
void Catalan() {
    f[0] = 1; f[1] = 1;
    for (int i=2; i<=100; ++i) {
        f[i] = f[i-1] * (4*i-2) / (i+1);
    }
}
