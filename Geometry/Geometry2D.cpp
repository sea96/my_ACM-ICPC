/*
 *点，直线，圆的定义，相交问题，面积问题，凸包等等
 */
const double EPS = 1e-10;
const double PI = acos(-1.0);

int sign(double x) {  //三态函数，减少精度问题
    return abs(x) < EPS ? 0 : x < 0 ? -1 : 1;
}
struct Point { //点的定义
    double x, y;
    Point(double x=0.0, double y=0.0) : x (x), y(y) {}
    Point operator + (const Point &rhs) const {  //向量加法
        return Point(x + rhs.x, y + rhs.y);
    }
    Point operator - (const Point &rhs) const {  //向量减法
        return Point(x - rhs.x, y - rhs.y);
    }
    Point operator * (double p) const {  //向量乘以标量
        return Point(x * p, y * p);
    }
    Point operator / (double p) const {  //向量除以标量
        return Point(x / p, y / p);
    }
    bool operator < (const Point &rhs) const {  //点的坐标排序
        return x < rhs.x || (x == rhs.x && y < rhs.y);
    }
    bool operator == (const Point &rhs) const {  //判断同一个点
        return sign(x - rhs.x) == 0 && sign(y - rhs.y) == 0;
    }
    void read() {
        scanf("%lf%lf", &x, &y);
    }
    void print() {
        printf("(%.2f, %.2f)\n", x, y);
    }
};
typedef Point Vector;  //向量的定义

double dot(Vector A, Vector B) {  //向量点积：向量长度的乘积再乘上夹角的余弦值，可判断两向量的角度
    return A.x * B.x + A.y * B.y;
}
double cross(Vector A, Vector B) {  //向量叉积，两向量组成的三角形的有向面积的两倍，可判断两向量的方向
    return A.x * B.y - A.y * B.x;
}
double polar_angle(Vector A) {  //向量极角
    return atan2(A.y, A.x);
}
double length(Vector A) {  //向量长度，点积
    return sqrt(dot(A, A));
}
double angle(Vector A, Vector B) {  //向量转角，逆时针，点积
    return acos(dot(A, B)/length(A)/length(B));
}
Vector rotate(Vector A, double rad) {  //向量旋转，逆时针
    return Vector(A.x*cos(rad)-A.y*sin(rad), A.x*sin(rad)+A.y*cos(rad));
}
Vector nomal(Vector A) {  //向量的单位法向量，确保A不是零向量
    double len = length(A);
    return Vector(-A.y/len, A.x/len);
}
Point line_line_inter(Point p, Vector V, Point q, Vector W) {  //两直线交点，参数方程
    Vector U = p - q;
    double t = cross(W, U) / cross(V, W);  //两直线有唯一交点的前提：cross(V, W)非0
    return p + V * t;
}
double point_to_line(Point p, Point a, Point b) {  //点到直线的距离，两点式
    Vector V1 = b - a, V2 = p - a;
    return fabs(cross(V1, V2)) / length(V1);  //点p到直线ab的距离
}
Point point_line_proj(Point p, Point a, Point b) {  //点在直线上的投影，两点式
    Vector V = b - a;
    return a + V * (dot (V, p-a)/dot(V, V));  //点p到直线ab的投影
}
bool can_seg_seg_inter(Point a1, Point a2, Point b1, Point b2) {  //判断线段规范相交，两点式
    double c1 = cross(a2-a1, b1-a1), c2 = cross(a2-a1, b2-a1),
           c3 = cross(b2-b1, a1-b1), c4 = cross(b2-b1, a2-b1);  //线段a1a2是否与线段b1b2规范相交
    return sign(c1)*sign(c2) < 0 && sign(c3)*sign(c4) < 0;  //若允许端点处相交，on_seg函数另外判断
}
bool on_seg(Point p, Point a, Point b) {  //判断点在线段上（不包含端点），两点式
    return sign(cross(a-p, b-p)) == 0 && sign(dot(a-p, b-p)) < 0;  //点p是否在线段ab上
}
double point_to_seg(Point p, Point a, Point b) {  //点到线段的距离，两点式
    if (a == b) return length(p - a);
    Vector V1 = b - a, V2 = p - a, V3 = p - b;  //点p到线段ab的距离
    if (sign(dot(V1, V2)) < 0) return length(V2);  //|pa|
    else if (sign(dot(V1, V3)) > 0) return length(V3);  //|pb|
    else return fabs(cross(V1, V2)) / length(V1);
}
bool can_line_seg_inter(Point a1, Point a2, Point b1, Point b2) {  //判断直线与线段相交（不包含端点），两点式
    double c1 = cross(a2-a1, b1-a1), c2 = cross(a2-a1, b2-a1);
    return sign(c1*c2) < 0;  //直线a1a2是否与线段b1b2相交，若允许线段端点处相交，改为<=
}
double area_triangle(Point a, Point b, Point c) {  //三角形面积，叉积
    return fabs(cross(b-a, c-a)) / 2.0;
}

typedef vector<Point> Polygon;  //多边形的定义。使用函数的前提：多边形poly排好序

//多边形的有向面积，叉积
double area_poly(Polygon &poly) {
    double ret = 0;
    for (int i=1; i<poly.size ()-1; ++i)  {
        ret += cross(poly[i]-poly[0], poly[i+1]-poly[0]);
    }
    return ret / 2.0;
}
//用有向直线a->b切割多边形poly，返回左侧点集
Polygon cut_polygon(Polygon &poly, Point a, Point b) {
    Polygon npoly;
    int n = poly.size();
    for (int i=0; i<n; ++i) {
        Point &c = poly[i];
        Point &d = poly[(i+1)%n];
        if (sign(cross(b-a, c-a)) >= 0) npoly.push_back(c);
        if (sign(cross(b-a, c-d)) != 0) {
            Point ip = line_line_inter(a, b-a, c, d-c);
            if (on_seg(ip, c, d)) npoly.push_back(ip);
        }
    }
    return npoly;  //如果退化，可能返回单点或线段
}
/*
 *判断点在多边形内，多边形可以为凹多边形甚至可以自交
 *凸多边形只需用叉积判断是否所有边在左边即可
 */
int point_in_polygon(Point p, Polygon &poly) {
    int wn = 0;
    int n = poly.size();
    for (int i=0; i<n; ++i) {
        if (on_seg(p, poly[i], poly[(i+1)%n]) || p == poly[i] || p == poly[(i+1)%n]) return -1;  //在多边形上
        int k = sign(cross(poly[(i+1)%n]-poly[i], p-poly[i]));
        int d1 = sign(poly[i].y-p.y);
        int d2 = sign(poly[(i+1)%n].y-p.y);
        if (k > 0 && d1 <= 0 && d2 > 0) wn++;
        if (k < 0 && d1 > 0 && d2 <= 0) wn--;
    }
    if (wn != 0) return 1;  //内部
    else return 0;  //外部
}
/*
 *凸包，基于水平序的Andrew算法。输入点的集合，返回凸包点的集合。
 *凸包边上无点：<=；凸包边上有点：<
 *时间复杂度O(NlogN)
 */
Polygon convex_hull(Polygon P) {
    sort(P.begin(), P.end());  //排序
    P.erase(unique(P.begin(), P.end()), P.end());  //删除重复点
    int n = P.size(), k = 0;
    Polygon Q(n*2);
    for (int i=0; i<n; ++i) {
        while (k > 1 && cross(Q[k-1]-Q[k-2], P[i]-Q[k-2]) <= 0) k--;
        Q[k++] = P[i];
    }
    for (int t=k, i=n-2; i>=0; --i) {
        while (k > t && cross(Q[k-1]-Q[k-2], P[i]-Q[k-2]) <= 0) k--;
        Q[k++] = P[i];
    }
    Q.resize(k-1);
    return Q;
}

struct Line {  //有向直线的定义，左边就是对应的半平面
    Point p;  //直线上任意一点
    Vector v;  //方向向量
    double ang;  //极角，即从x正方向旋转到向量v所需的角（弧度）
    Line() {}
    Line(Point p, Vector v) : p(p), v(v) {
        ang = polar_angle(v);
    }
    bool operator < (const Line &rhs) const {
        return ang < rhs.ang;
    }
    Point point(double a) {
        return p + v * a;
    }
};

/*
 *半平面交，输入直线集，返回凸多边形的点集，如果不存在返回的是空集
 *时间复杂度O(NlogN)
 */
bool point_on_left(Point p, Line L) {  //判断点在直线左侧
    return cross(L.v, p-L.p) > 0;
}
Polygon half_plane_inter(vector<Line> L) {
    sort(L.begin(), L.end());  //按极角排序
    int first, last;  //双端队列的第一个元素和最后一个元素的下标
    int n = L.size();
    Point *p = new Point[n];  //p[i]为q[i]与q[i+1]的交点
    Line *q = new Line[n];  //双端队列，保存半平面
    q[first=last=0] = L[0];  //初始化只有一个半平面L[0]
    for (int i=1; i<n; ++i) {
        while (first < last && !point_on_left(p[last-1], L[i])) last--;
        while (first < last && !point_on_left(p[first], L[i])) first++;
        q[++last] = L[i];
        if (sign(cross(q[last].v, q[last-1].v)) == 0) {
            last--;
            if (point_on_left(L[i].p, q[last])) q[last] = L[i];
        }
        if (first < last) p[last-1] = line_line_inter(q[last-1].p, q[last-1].v, q[last].p, q[last].v);
    }
    while (first < last && !point_on_left(p[last-1], q[first])) last--;
    Polygon poly;
    //删去无用平面
    if (last - first <= 1)  return poly;  //空集
    p[last] = line_line_inter(q[last].p, q[last].v, q[first].p, q[first].v);  //首尾两个半平面的交点
    for (int i=first; i<=last; ++i) poly.push_back(p[i]);  //复制结果
    return poly;
}

struct Circle {  //圆的定义
    Point c;  //圆心
    double r;  //半径
    Circle() {}
    Circle(Point c, double r) : c(c), r(r) {}
    Point point(double a) {  //圆上的一点
        return Point(c.x+cos(a)*r, c.y+sin(a)*r);
    }
};

//直线与圆相交求交点，返回交点个数，交点保存在P中
int line_cir_inter(Line L, Circle C, double &t1, double &t2, vector<Point> &P) {
    double a = L.v.x, b = L.p.x-C.c.x, c = L.v.y, d = L.p.y-C.c.y;
    double e = a*a+c*c, f = 2.0*(a*b+c*d), g = b*b+d*d-C.r*C.r;
    double delta = f*f-4.0*e*g;  //判别式
    //以上参数的具体含义参见《训练指南》P264
    //交点p1=L.p+L.v*t1，交点p2=L.p+L.v*t2
    if (sign(delta) < 0) return 0;  //相离
    if (sign(delta) == 0) {  //相切
        t1 = t2 = -f / (2.0*e);
        P.push_back(L.point(t1));
        return 1;
    }
    //相交
    t1 = (-f-sqrt(delta)) / (2.0*e);
    t2 = (-f+sqrt(delta)) / (2.0*e);
    if (t1 > t2) swap (t1, t2);
    if (sign(t1) > 0 && sign(t1-1.0) < 0) P.push_back(L.point(t1));
    if (sign(t2) > 0 && sign(t2-1.0) < 0) P.push_back(L.point(t2));
    return 2;
}

//两圆相交求交点，返回交点个数，交点保存在P中
int cir_cir_inter_point(Circle C1, Circle C2, vector<Point> &P) {
    double d = length(C1.c - C2.c);
    if (sign(d) == 0) {  //同心圆
        if (sign(C1.r-C2.r) == 0) return -1;  //两圆重合
        else return 0;
    }
    if (sign(C1.r+C2.r-d) < 0) return 0; //外离
    if (sign(fabs(C1.r-C2.r)-d) < 0) return 0;  //内含
    double a = polar_angle(C2.c - C1.c);
    double da = acos((C1.r*C1.r+d*d-C2.r*C2.r) / (2*C1.r*d));  //C1C2到C1P1的角？
    Point p1 = C1.point(a - da), p2 = C2.point(a + da);
    P.push_back(p1);
    if (p1 == p2) return 1;
    else P.push_back(p2);
    return 2;
}
//余弦定理，abc三条边，返回c边对应的角C的cos值
double squ(double x) {  //平方
    return x * x;
}
double distance(Point a, Point b) {  //两点的距离
    return sqrt(squ(a.x-b.x)+squ(a.y-b.y));
}
double cosine(double a, double b, double c) {  //输入三角形的三条边，角C的余弦值
    return (squ(a)+squ(b)-squ(c)) / (2*a*b);
}
//两圆相交求相交面积，返回面积值
double cir_cir_inter_area(Circle C1, Circle C2) {
    double dist = distance(C1.c, C2.c);
    if (C1.r + C2.r < dist) return 0.0;  //相离
    else if (sign(abs(C1.r-C2.r)-dist) >= 0) {  //内含
        int r = C1.r < C2.r ? C1.r : C2.r;
        return PI*r*r;
    } else {  //相交
        double ang1 = 2.0 * acos(cosine(C1.r, dist, C2.r));
        double ang2 = 2.0 * acos(cosine(C2.r, dist, C1.r));
        double ret = squ(C1.r)*ang1/2 + squ(C2.r)*ang2/2 - squ(C1.r)*sin(ang1)/2 - squ(C2.r)*sin(ang2)/2;  //面积容斥
        return ret;
    }
}

//过点到圆的切线，返回切线条数，切线的向量保存在V中
int point_cir_tan(Point p, Circle C, Vector* V) {
    Vector u = C.c - p;
    double dis = length(u);
    if (dis < C.r) return 0;  //点在圆内
    else if (sign(dis-C.r) == 0) {  //点在圆上
        V[0] = rotate(u, PI/2);
        return 1;
    } else {  //点在圆外
        double ang = asin(C.r/dis);
        V[0] = rotate(u, -ang);
        V[1] = rotate(u, +ang);
        return 0;
    }
}
//两圆的公切线，返回公切线条数，a[i]和b[i]分别表示第i条切线在圆A和圆B的切点
int cir_cir_tan(Circle A, Circle B, Point* a, Point* b) {
    if (A.r < B.r) {
        swap(A, B);
        swap(a, b);
    }
    double d = dot(A.c-B.c, A.c-B.c);
    double rsub = A.r-B.r, rsum = A.r+B.r;
    if (sign(d-rsub) < 0) return 0; //内含，零条公切线
    double base = polar_angle(B.c - A.c);
    if (sign(d) == 0 && sign(A.r-B.r) == 0) return -1; //两圆重合，无数条公切线
    if (sign(d-rsub) == 0) {  //内切，一条公切线
        a[0] = A.point(base); b[0] = B.point(base);
        return 1;
    }
    //有外公切线
    int cnt = 0;
    double ang = acos(rsub/d);
    a[cnt] = A.point(base+ang); b[cnt] = B.point(base+ang); cnt++;
    a[cnt] = A.point(base-ang); b[cnt] = B.point(base-ang); cnt++;
    if (d == rsum) {  //一条内公切线
        a[cnt] = A.point (base); b[cnt] = B.point(base + PI); cnt++;
    } else if (sign (d-rsum) > 0) {  //两条内公切线
        double ang2 = acos(rsum / d);
        a[cnt] = A.point(base+ang2); b[cnt] = B.point(base+ang2+PI); cnt++;
        a[cnt] = A.point(base-ang2); b[cnt] = B.point(base-ang2+PI); cnt++;
    }
    return cnt;
}
/*
   多边形与圆的公共面积，上交红书模板
   调用fabs (cir_poly_area (ps))，ps为多边形的点集，
   需要用到line_cir_inter ()函数，圆心在原点(可平移)
 */
double sector_area(Point a, Point b, double r)    {     //三角剖分，求扇形面积
    double theta = polar_angle (a) - polar_angle (b);
    while (sign (theta) <= 0)   theta += 2 * PI;
    while (theta > 2 * PI)  theta -= 2 * PI;
    theta = min (theta, 2 * PI - theta);
    return r * r * theta / 2;
}
double cal(Point a, Point b, double r)    {
    double t1, t2;
    bool ina = sign (length (a) - r) < 0;
    bool inb = sign (length (b) - r) < 0;
    if (ina && inb) return fabs (cross (a, b)) / 2.0;
    vector<Point> p;
    int num = line_cir_inter (Line (a, b - a), Circle (Point (0, 0), r), t1, t2, p);
    if (ina)    return sector_area (b, p[0], r) + fabs (cross (a, p[0])) / 2.0;
    if (inb)    return sector_area (p[0], a, r) + fabs (cross (p[0], b)) / 2.0;
    if (num == 2)   return sector_area (a, p[0], r) + sector_area (p[1], b, r) + fabs (cross (p[0], p[1])) / 2.0;
    return sector_area (a, b, r);
}
double cir_poly_area(vector<Point> &ps, double r)  {
    double ret = 0;
    for (int i=0; i<ps.size ()-1; ++i) {        //多边形最后放入ps[0]起点
        int sgn = sign (cross (ps[i], ps[i+1]));
        if (sgn != 0)   {
            ret += sgn * cal (ps[i], ps[i+1], r);
        }
    }
    return ret;
}
/*
 *最近点对问题，分治法，先按照x坐标排序，求解(left, mid)和(mid+1, right)范围的最小值，
 *然后类似区间合并，分离mid左右的点也求最小值
 */
double min_dist(int left, int right) {
    if (left == right) return INF;
    else if (right - left == 1) return get_dist(point[left], point[right]);
    else {
        int mid = left + right >> 1;
        double ret = min(min_dist(left, mid), min_dist(mid + 1, right));
        if (ret == 0) return ret;
        int endy = 0;
        for (int i=mid; i>=left&&point[mid].x-point[i].x<=ret; --i) idy[endy++] = i;
        for (int i=mid+1; i<=right&&point[i].x-point[mid+1].x<=ret; ++i) idy[endy++] = i;
        sort (idy, idy+endy, cmp_y);  //return point[i].y < point[j].y;
        for (int i=0; i<endy; ++i) {
            for (int j=i+1; j<endy&&point[j].y-point[i].y<ret; ++j) {
                ret = min(ret, get_dist(point[i], point[j]));
            }
        }
        return ret;
    }
}
