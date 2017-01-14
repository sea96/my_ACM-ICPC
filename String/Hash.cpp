const int BASE = 233;
unsigned long long Hash[N], pw[N];
void calcHash() {
    pw[0] = 1;
    for (int i=1; i<len; ++i) pw[i] = pw[i-1] * BASE;
    Hash[len] = 0;
    for (int i=len-1; i>=0; --i)
        Hash[i] = Hash[i+1]*BASE+(s[i]-'a');
    //Hash[i, i+L) = H[i] - H[i+L]*pw[L];
}
