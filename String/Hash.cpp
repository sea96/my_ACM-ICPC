void get_hash(char *s) {
    int len = strlen(s);
    H[len] = 0;
    for (int i=len-1; i>=0; --i) {
        H[i] = H[i+1] * x + (s[i] - 'a');
    }
    xp[0] = 1;
    for (int i=1; i<len; ++i) {
        xp[i] = xp[i-1] * x;
    }
    //hash[i, i+L) = H[i] - H[i+L] * xp[L];
}
