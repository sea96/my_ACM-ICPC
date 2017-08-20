// 用画有限自动机理解
void getFail(char *p, int len) {
    int j = fail[0] = fail[1] = 0;
    for (int i=1; i<len; ++i) {
        while (j && p[i] != p[j]) j = fail[j];
        fail[i+1] = p[i] == p[j] ? ++j : 0;
    }
}

int KMP(char *T, char *P) {
    int len_t = strlen(T), len_p = strlen(P);
    getFail(P, len_p);
    for (int j=0, i=0; i<len_t; ++i) {
        while (j && P[j] != T[i]) j = fail[j];
        j += (P[j] == T[i]);
        if (j == lenp) return i - j + 1;
        //统计出现次数：ans++; j = 0（不重复）/fail[j]（可重复）
    }
    return -1;
}
