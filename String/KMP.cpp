void getFail(char *p, int len) {
    fail[0] = fail[1] = 0;
    for (int i=1; i<len; ++i) {
        int j = fail[i];
        while (j && p[i] != p[j]) j = fail[j];
        fail[i+1] = p[i] == p[j] ? j+1 : 0;
    }
}

int KMP(char *T, char *P) {
    int lent = strlen(T);
    int lenp = strlen(P);
    get_fail(P, lenp);
    int i = 0, j = 0;
    while (i < lent) {
        while (j != -1 && T[i] != P[j]) j = fail[j];
        i++; j++;
        if (j == lenp) return i - j + 1;
        //统计出现次数：ans++; j = 0（不重复）/fail[j]（可重复）
    }
    return -1;
}
