/*
    *后缀数组，倍增算法实现，复杂度O(nlogn)
    *sa[i]: 第i小的后缀是在字符串位置，即后缀sa[i]
    *rank[i]: 后缀i在sa数组下标，即第rank[i]小
    *height[i]: LCP (suffix (sa[i-1], sa[i]))
*/
struct Suffix_Array {
    int n, len, s[N];
    int sa[N], rank[N], height[N];
    int tmp_one[N], tmp_two[N], c[N];
    
    void init_str(char *str);
    void build_sa(int m = 128);
    void get_height();

	void print();
	
	void RMQ_init();
	int RMQ_query(int l, int r);
};

//RMQ_query (rank[i], rank[j]);
int Suffix_Array::RMQ_query(int l, int r) {
    if (l > r) {
        std::swap (l, r);
    }
    l++;
    int k = 0; while (1<<(k+1) <= r - l + 1) k++;
    return std::min (dp[l][k], dp[r-(1<<k)+1][k]);
}

void Suffix_Array::print(char *str) {
    puts ("/*Suffix*/");
    for (int i=0; i<n; ++i) {
        printf ("%s\n", str+sa[i]);
    }
}

void Suffix_Array::RMQ_init() {
    for (int i=0; i<n; ++i) {
        dp[i][0] = height[i];
    }
    for (int j=1; (1<<j)<=n; ++j) {
        for (int i=0; i+(1<<j)-1<n; ++i) {
            dp[i][j] = std::min (dp[i][j-1], dp[i+(1<<(j-1))][j-1]);
        }
    }
}

void Suffix_Array::init_str(char *str) {
    n = 0;
    len = strlen (str);
    for (int i=0; i<len; ++i) {
        s[n++] = str[i] - 'a' + 1;
    }
    s[n++] = 0;  //n = strlen (str) + 1
}

void Suffix_Array::get_height() {
    for (int i=0; i<n; ++i) rank[sa[i]] = i;
    int k = height[0] = 0;
    for (int i=0; i<n-1; ++i) {
        if (k) k--;
        int j = sa[rank[i]-1];
        while (s[i+k] == s[j+k]) k++;
        height[rank[i]] = k;
    }
}

//m = max (r[i]) + 1，一般字符128足够了
void Suffix_Array::build_sa(int m) {
    int i, j, p, *x = tmp_one, *y = tmp_two;
    for (i=0; i<m; ++i) c[i] = 0;
    for (i=0; i<n; ++i) c[x[i]=s[i]]++;
    for (i=1; i<m; ++i) c[i] += c[i-1];
    for (i=n-1; i>=0; --i) sa[--c[x[i]]] = i;
    for (j=1, p=1; p<n; j<<=1, m=p) {
        for (p=0, i=n-j; i<n; ++i) y[p++] = i;
        for (i=0; i<n; ++i) if (sa[i] >= j) y[p++] = sa[i] - j;
        for (i=0; i<m; ++i) c[i] = 0;
        for (i=0; i<n; ++i) c[x[y[i]]]++;
        for (i=1; i<m; ++i) c[i] += c[i-1];
        for (i=n-1; i>=0; --i) sa[--c[x[y[i]]]] = y[i];
        std::swap (x, y);
        for (p=1, x[sa[0]]=0, i=1; i<n; ++i) {
            x[sa[i]] = (y[sa[i-1]] == y[sa[i]] && y[sa[i-1]+j] == y[sa[i]+j] ? p - 1 : p++);
        }
    }
}
