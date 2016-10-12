/*
    T[]是文本串，P[]是模式串，lent，lenp是各自的长度
    返回第一个匹配完成在文本串上的位置，不成功返回-1
*/
void get_fail(char *P, int lenp) {
    int i = 0, j = -1;
    fail[0] = -1;
    while (i < lenp) {
        if (j == -1 || P[j] == P[i]) {
            i++; j++; fail[i] = j;
        }
        else j = fail[j];
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
