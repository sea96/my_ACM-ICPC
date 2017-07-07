//http://www.cnblogs.com/yefeng1627/archive/2012/12/24/2830979.html
char S[N], T[N];
int lcp_t[N];  // T的后缀i和T的最长公共前缀长度
int lcp_s[N];  // S的后缀i和T的最长公共前缀长度
int len_s, len_t;

void getLCP(char *T) {
    lcp_t[0] = len_t;
    int i, a = 1;  // [0,i)里能匹配最远位置的点
    for(i = 0; i<len_t-1 && T[i]==T[i+1]; ++i);
    lcp_t[1] = i;
    for (i=2; i<len_t; ++i) {
        int p = a+lcp_t[a]-1, L = lcp_t[i-a];  // i在另一个T的对应点是i-a
        if (i+L-1 >= p) {  // 包括等于的情况
            int j = max(p-i+1, 0);  // 可能发生 p-i+1<0 的情况
            while (i+j < len_t && T[i+j]==T[j]) j++;  // 从 p+1 和 p-i+1 开始往后继续匹配
            lcp_t[i] = j; a = i;
        } else lcp_t[i] = L;
    }
}

void extendedKMP(char *S, char *T) {
    getLCP(T);
    int i = 0, a = 0;
    while (i<len_s && i<len_t && S[i] == T[i]) i++;
    lcp_s[0] = i;
    for (i=1; i<len_s; ++i) {
        int p = a+lcp_s[a]-1, L = lcp_t[i-a];
        if(i+L-1 >= p) {
            int j = max(p-i+1, 0);
            while (i+j < len_s && j<len_t && S[i+j] == T[j]) j++;
            lcp_s[i] = j; a = i;
        } else lcp_s[i] = L;
    }
}
