/*
 *Treap名次树，优先级随机，插入，删除，查找期望时间复杂度均为O(logn)
 *主要功能：求第k大元素，求x的名次(比x小的个数+1)
 */
#include <bits/stdc++.h>

//Treap
struct Node {
    Node* ch[2];
    int v, r, s, c;
    Node(int v) : v(v) {
        ch[0] = ch[1] = NULL;
        r = rand();
        s = c = 1;
    }
    bool operator < (const Node& rhs) const {
        return r < rhs.r;
    }
    int cmp(int x) const {
        if (x == v) return -1;
        return x < v ? 0 : 1;
    }
    void maintain() {
        s = c;
        if (ch[0]) s += ch[0]->s;
        if (ch[1]) s += ch[1]->s;
    }
};

Node* root;

void rotate(Node* &o, int d) {
    Node* k = o->ch[d^1]; o->ch[d^1] = k->ch[d]; k->ch[d] = o;
    o->maintain(); k->maintain(); o = k;
}

void insert(Node* &o, int x) {
    if (!o) o = new Node(x);
    else if (o->v == x) o->c += 1;
    else {
        int d = o->v < x ? 1 : 0;
        insert(o->ch[d], x);
        if (o->ch[d] > o) rotate(o, d^1);
    }
    o->maintain();
}

void remove_node(Node* &o, int x) {
    int d = o->cmp(x);
    if (d == -1) {
        if (o->c > 1) o->c -= 1;
        else {
            Node* u = o;
            if (o->ch[0] && o->ch[1]) {
                int d2 = o->ch[0] > o->ch[1] ? 1 : 0;
                rotate(o, d2);
                remove_node(o->ch[d2], x);
            } else {
                if (o->ch[0]) o = o->ch[0];
                else o = o->ch[1];
                delete u;
            }
        }
    } else remove_node(o->ch[d], x);
    if (o) o->maintain();
}

int kth(Node* &o, int k) {
    if (!o || k <= 0 || k > o->s) return -1;
    int lsum = o->ch[0]==NULL ? 0 : o->ch[0]->s;
    if (k <= lsum) return kth(o->ch[0], k);
    else if (k > lsum + o->c) return kth(o->ch[1], k-lsum-o->c);
    else return o->v;
}

int rank(Node* &o, int x) {
    if (!o) return 0;
    int lsum = o->ch[0]==NULL ? 0 : o->ch[0]->s;
    if (x == o->v) return lsum + 1;
    if (x < o->v) return rank(o->ch[0], x);
    else return lsum + o->c + rank(o->ch[1], x);
}

void query_pro(Node* &o, int x, int &ans) {
    if (!o) return ;
    if (o->v < x) {
        ans = o->v;
        query_pro(o->ch[1], x, ans);
    } else query_pro(o->ch[0], x, ans);
}

void query_sub(Node* &o, int x, int &ans) {
    if (!o) return ;
    if (o->v > x) {
        ans = o->v;
        query_sub(o->ch[0], x, ans);
    } else query_sub(o->ch[1], x, ans);
}

void remove_tree(Node* &o)  {
    if (!o) return ;
    remove_tree(o->ch[0]);
    remove_tree(o->ch[1]);
    delete o;
    o = NULL;
}
int main() {
    int n;
    while (scanf("%d", &n) == 1) {
        root = NULL;
        int opt, x, ans;
        for (int i=1; i<=n; ++i) {
            scanf("%d%d", &opt, &x);
            if (opt == 1) insert(root, x);
            else if (opt == 2) remove_node(root, x);
            else if (opt == 3) printf("%d\n", rank(root, x));
            else if (opt == 4) printf("%d\n", kth(root, x));
            else if (opt == 5) {
                ans = -1;
                query_pro(root, x, ans);
                printf("%d\n", ans);
            }
            else if (opt == 6) {
                ans = -1;
                query_sub(root, x, ans);
                printf("%d\n", ans);
            }
        }
        remove_tree(root);
    }
    return 0;
}
//
struct Node {
    Node* ch[2]; //左右孩子
    int v, r, s; //值，优先级，子树数量
    Node(int v) : v (v) {
        ch[0] = ch[1] = NULL;   r = rand ();    s = 1;
    }
    bool operator < (const Node& rhs) const {
        return r < rhs.r;
    }
    int cmp(int x) const {
        if (x == v) return -1;
        return x < v ? 0 : 1;
    }
    void maintain() {
        s = 1;
        if (ch[0])  s += ch[0]->s;
        if (ch[1])  s += ch[1]->s;
    }
};
//旋转，d=0表示左旋，d=1表示右旋
void rotate(Node* &o, int d) {
    Node* k = o->ch[d^1];   o->ch[d^1] = k->ch[d];  k->ch[d] = o;
    o->maintain (); k->maintain (); o = k;
}
/*
 *插入，如果x比o值小，插入左子树，否则插入右子树。
 *插入完毕后，如果优先级破坏，则往反方向旋转
 */
void insert(Node* &o, int x) {
    if (!o) o = new Node (x);
    else {
        int d = (x < o->v ? 0 : 1);
        insert (o->ch[d], x);
        if (o->ch[d]->r > o->r) rotate (o, d ^ 1);
    }
    o->maintain ();
}
/*
 *删除，先找到待删除结点，如果只有一棵子树，替换即可
 *若有两棵子树，则将子树中优先级高的旋转到根，然后递归另一棵子树去删除
 */
void remove_node(Node* &o, int x) {
    int d = o->cmp (x);
    if (d == -1) {
        Node* u = o;
        if (o->ch[0] && o->ch[1]) {
            int d2 = (o->ch[0]->r > o->ch[1]->r ? 1 : 0);
            rotate (o, d2); remove_node (o->ch[d2], x);
        }
        else {
            if (o->ch[0])   o = o->ch[0];
            else    o = o->ch[1];
            delete u;
        }
    }
    else    remove_node (o->ch[d], x);
    if (o)  o->maintain ();
}
void remove_tree(Node* &o) {
    if (o->ch[0])   remove_tree (o->ch[0]);
    if (o->ch[1])   remove_tree (o->ch[1]);
    delete o;
    o = NULL;
}
//第k大的值，注意是第k大！
int kth(Node* &o, int k) {
    if (!o || k <= 0 || k > o->s)   return 0;
    int sum = (o->ch[1] == NULL ? 0 : o->ch[1]->s);
    if (k == sum + 1)   return o->v;
    else if (k <= sum)  return kth (o->ch[1], k);
    else    return kth (o->ch[0], k - sum - 1);
}
void merge(Node* &src, Node* &dest) {
    if (src->ch[0]) merge (src->ch[0], dest);
    if (src->ch[1]) merge (src->ch[1], dest);
    insert (dest, src->v);
    delete src;
    src = NULL;
}
