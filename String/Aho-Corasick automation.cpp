//AC自动机
struct AC {
    static const int NODE = 2000;
    static const int SIZE = 26;
    //ch[][]树上结点，fail[]指针，last[]沿着失配边往回走遇到的下一个结点编号
    int ch[NODE][SIZE], fail[NODE], last[NODE];
    int end[NODE];
    int sz;
    
    void clear() {
        memset (ch[0], 0, sizeof (ch[0]));
        end[0] = 0;
        sz = 1;
    }
    int idx(char ch) {
        return ch - 'a';
    }
    void insert(char *str, int v) {
        int u = 0;
        for (int c, i=0; str[i]; ++i) {
            c = idx (str[i]);
            if (!ch[u][c]) {
                memset (ch[sz], 0, sizeof (ch[sz]));
                end[sz] = 0;
                ch[u][c] = sz++;
            }
            u = ch[u][c];
        }
        end[u] = v;
    }
    void build() {
        fail[0] = 0;
        std::queue<int> que;
        for (int c=0; c<SIZE; ++c) {
            int u = ch[0][c];
            if (u) {
                fail[u] = 0;
                last[u] = 0;
                que.push (u);
            }
        }
        while (!que.empty ()) {
            int r = que.front (); que.pop ();
            for (int c=0; c<SIZE; ++c) {
                int u = ch[r][c];
                if (!u) {
                    ch[r][c] = ch[fail[r]][c];
                } else {
                    fail[u] = ch[fail[r]][c];
                    last[u] = end[fail[u]] ? fail[u] : last[fail[u]];
                    que.push (u);
                }
            }
        }
    }
};
