/*
    T[]是文本串，P[]是模式串，lent，lenp是各自的长度
    返回第一个匹配完成在文本串上的位置，不成功返回-1
*/
void get_fail(char *P) {
    int i = 0, j = -1;
    fail[0] = -1;
    while (i < lenP) {
        if (j == -1 || P[j] == P[i]) {
            i++; j++; fail[i] = j;
        }
        else j = fail[j];
    }
}

void KMP(char *T, char *P) {
    lenT = strlen (T);
    lenP = strlen (P);
    get_fail (P);
    int i = 0, j = 0;
    while (i < lenT) {
        while (j != -1 && T[i] != P[j]) j = fail[j];
        i++; j++;
        if (j == lenP) return i - j + 1;
    }
}
