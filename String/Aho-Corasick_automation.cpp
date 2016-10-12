struct AC {
    static const int NODE = N*L;  //整个Trie的节点数=单词个数*单词长度
    static const int SIZE = 26;  //单词字符集合
    //节点，失配指针（“备胎”），下一个结点，是否为结点（结束的点）
    int ch[NODE][SIZE], fail[NODE], last[NODE], end[NODE];
    int sz;  //当前Trie的节点数

    void init();
    int idx(char c);
    void insert(char *s);
    void get_fail();
}ac;

void AC::init() {
    memset(ch[0], 0, sizeof (ch[0]));
    end[0] = 0;
    sz = 1;
}

int AC::idx(char c) {
    return c - 'a';
}

void AC::insert(char *s) {
    int u = 0;
    for (int c, i=0; s[i]; ++i) {
        c = idx(s[i]);
        if (!ch[u][c]) {
            memset(ch[sz], 0, sizeof (ch[sz]));
            end[sz] = 0;
            ch[u][c] = sz++;
        }
        u = ch[u][c];
    }
    end[u] = 1;  //结束标记
}

void AC::get_fail() {
    fail[0] = last[0] = 0;
    queue<int> que;
    for (int c=0; c<SIZE; ++c) {
        int u = ch[0][c];
        if (u) {
            fail[u] = 0;
            last[u] = 0;
            que.push (u);
        }
    }
    while (!que.empty ()) {
        int r = que.front ();
        que.pop ();
        for (int c=0; c<SIZE; ++c) {
            int &u = ch[r][c];
            if (!u) {
                //一视同仁，把不存在的边补上，但并不保证转移有效
                u = ch[fail[r]][c];
            } else {
                int v = fail[r];
                while (v && !ch[v][c]) v = fail[v];
                fail[u] = ch[v][c];  //找到“备胎”的位置（也不保证有）
                end[u] |= end[fail[u]];  //提高效率，是否存在匹配的模板串
                last[u] = end[fail[u]] ? fail[u] : last[fail[u]];  //找“备胎”中下一个结点
                que.push (u);
            }
        }
    }
    //事实上，调用函数后，再也不需要失配指针了
}
