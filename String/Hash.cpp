const int BASE = 233;
unsigned long long _hash[N], pw[N];
void calcHash() {
    pw[0] = 1;
    for (int i=1; i<len; ++i) pw[i] = pw[i-1] * BASE;
    _hash[len] = 0;
    for (int i=len-1; i>=0; --i)
        _hash[i] = _hash[i+1]*BASE+(s[i]-'a');
    //_hash[i, i+L) = _hash[i] - _hash[i+L]*pw[L];
}
