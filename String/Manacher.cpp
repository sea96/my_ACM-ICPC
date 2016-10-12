/*
 *最长回文串算法，O(N)
 */
int Manacher(char *s) {
    int n = strlen (s);
    //"abcde" -> "$#a#b#c#d#e#"
    for (int i=n; i>=0; --i) {
        s[i*2+2] = s[i];
        s[i*2+1] = '#';
    }
    s[0] = '$';
    n = n * 2 + 2;
    int id = 0; p[0] = 1;
    for (int i=2; i<n; ++i) {
        if (id + p[id] > i) p[i] = min(p[2*id-i], id + p[id] - i);
        else p[i] = 1;
        while (s[i-p[i]] == s[i+p[i]]) p[i]++;
        if (id + p[id] < i + p[i]) id = i;
    }
    int ret = 0;
    for (int i=0; i<n-1; ++i) {
        ret = max(ret, p[i] - 1);
    }
    return ret;
}
