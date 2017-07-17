typedef unsigned long long ull;
ull pw[N];
ull H[N];
ull seed = 233;
// ull容易被卡(CF822E)，但换MOD速度慢
void initHash(char *s, int len) {
    H[len] = 0;
    for (int i=len-1; i>=0; --i) {
        H[i] = H[i+1] * seed + s[i];
    }
    pw[0] = 1;
    for (int i=1; i<len; ++i) pw[i] = pw[i-1] * seed;
    //hash[i, i+L) = H[i] - H[i+L] * pw[L];
}
