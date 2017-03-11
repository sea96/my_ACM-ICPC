const int N = 1e5 + 5;
struct Node {
    Node *ch[2], *fa;  //儿子，父亲
    int v, s, flip;  //节点的值，子树个数，反转标记
    Node() { s = 0; }  //给null初始化
    int d() {
        return fa->ch[1]==this;  //判断自己是哪个儿子
    }
    void setc(Node* c, int d) {
        ch[d] = c; c->fa = this;  //指定儿子
    }
    void up() {
        s = ch[0]->s + ch[1]->s + 1;  //向上更新
    }
    void down() {
        if (flip) {  //向下更新
            flip ^= 1;
            swap(ch[0], ch[1]);
            ch[0]->flip ^= 1;
            ch[1]->flip ^= 1;
        }
    }
}Tnull, *null = &Tnull;
Node pool[N], *node = pool;  //多组测试时，指针每次要移到0

Node* new_node(int v) {
    Node *o = node++;
    o->ch[0] = o->ch[1] = o->fa = null;
    o->v = v; o->s = 1; o->flip = 0;
    return o;
}

//将节点o向上旋转
void rotate(Node* o) {
    Node* p = o->fa;
    int d = o->d();
    p->fa->setc(o, p->d());
    p->setc(o->ch[d^1], d);
    o->setc(p, d^1);
    p->up();
}

//将节点o伸展到p的儿子的位置
void splay(Node* o, Node *p = null) {
    while (o->fa != p) {
        if (o->fa->fa == p) rotate(o);
        else {
            //三点共线：先旋转父亲，再旋转o；不共线：o旋转两次
            if (o->d() == o->fa->d()) rotate(o->fa), rotate(o);
            else rotate(o), rotate(o);
        }
    }
    o->up();
}

//找到以o为根的树上第k小
Node* kth(Node* o, int k) {
    while (1) {
        o->down();  //向下走之前将懒惰标记向下传
        int c = 1 + o->ch[0]->s;
        if (k == c) return o;
        if (k > c) { k -= c; o = o->ch[1]; }
        else o = o->ch[0];
    }
    return null;
}

//将前k小放到左子树，其它的放到右子树
void split(Node* o, int k, Node* &left, Node* &right) {
    if (k == 0) { left = null; right = o; }
    else if (k == o->s) { left = o; right = null; }
    else {
        //先找到分裂节点，伸展到根再分裂
        right = kth(o, k+1);
        splay(right);
        left = right->ch[0];
        left->fa = null;
        right->ch[0] = right->fa = null;
        right->up();
    }
}

void merge(Node* &o, Node* left, Node* right) {
    if (left == null) o = right;
    else if (right == null) o = left;
    else {
        //指定左子树最后一个节点为根节点
        o = kth(left, left->s);
        splay(o);
        o->setc(right, 1);
        o->up();
    }
}
