/*************************************************/
计算几何
/*
	点，直线，圆的定义，相交问题，面积问题，凸包 等等
*/
const double EPS = 1e-10;
const double PI = acos (-1.0);
int dcmp(double x) {  //三态函数，减少精度问题
    if (fabs (x) < EPS) return 0;
    else    return x < 0 ? -1 : 1;
}
struct Point {  //点的定义
    double x, y;
    Point (double _x=0.0, double _y=0.0) : x (_x), y (_y) {}
    Point operator + (const Point &r) const {  //向量加法
        return Point (x + r.x, y + r.y);
    }
    Point operator - (const Point &r) const {  //向量减法
        return Point (x - r.x, y - r.y);
    }
    Point operator * (double p) const {  //向量乘以标量
        return Point (x * p, y * p);
    }
    Point operator / (double p) const {  //向量除以标量
        return Point (x / p, y / p);
    }
    bool operator < (const Point &r) const {  //点的坐标排序
        return x < r.x || (x == r.x && y < r.y);
    }
    bool operator == (const Point &r) const {  //判断同一个点
        return dcmp (x - r.x) == 0 && dcmp (y - r.y) == 0;
    }
    void read() {
        double _x, _y; scanf ("%lf%lf", &_x, &_y);
        x = _x; y = _y;
    }
    void print() {
        printf ("%.2f %.2f\n", x, y);
    }
};
typedef Point Vector;  //向量的定义

double dot(Vector A, Vector B)  {       //向量点积
    return A.x * B.x + A.y * B.y;
}
double cross(Vector A, Vector B)    {       //向量叉积
    return A.x * B.y - A.y * B.x;
}
double polar_angle(Vector A)  {     //向量极角
    return atan2 (A.y, A.x);
}
double length(Vector A) {       //向量长度，点积
    return sqrt (dot (A, A));
}
double angle(Vector A, Vector B)    {       //向量转角，逆时针，点积
    return acos (dot (A, B) / length (A) / length (B));
}
Vector rotate(Vector A, double rad) {       //向量旋转，逆时针
    return Vector (A.x * cos (rad) - A.y * sin (rad), A.x * sin (rad) + A.y * cos (rad));
}
Vector nomal(Vector A)  {       //向量的单位法向量
    double len = length (A);
    return Vector (-A.y / len, A.x / len);
}
Point line_line_inter(Point p, Vector V, Point q, Vector W)    {        //两直线交点，参数方程
    Vector U = p - q;
    double t = cross (W, U) / cross (V, W);
    return p + V * t;
}
double point_to_line(Point p, Point a, Point b)   {       //点到直线的距离，两点式
    Vector V1 = b - a, V2 = p - a;
    return fabs (cross (V1, V2)) / length (V1);
}
double point_to_seg(Point p, Point a, Point b)    {       //点到线段的距离，两点式
    if (a == b) return length (p - a);
    Vector V1 = b - a, V2 = p - a, V3 = p - b;
    if (dcmp (dot (V1, V2)) < 0)    return length (V2);
    else if (dcmp (dot (V1, V3)) > 0)   return length (V3);
    else    return fabs (cross (V1, V2)) / length (V1);
}
Point point_line_proj(Point p, Point a, Point b)   {     //点在直线上的投影，两点式
    Vector V = b - a;
    return a + V * (dot (V, p - a) / dot (V, V));
}
bool can_seg_seg_inter(Point a1, Point a2, Point b1, Point b2)  {       //判断线段相交，两点式
    double c1 = cross (a2 - a1, b1 - a1), c2 = cross (a2 - a1, b2 - a1),
           c3 = cross (b2 - b1, a1 - b1), c4 = cross (b2 - b1, a2 - b1);
    return dcmp (c1) * dcmp (c2) < 0 && dcmp (c3) * dcmp (c4) < 0;
}
bool can_line_seg_inter(Point a1, Point a2, Point b1, Point b2)    {        //判断直线与线段相交，两点式
    double c1 = cross (a2 - a1, b1 - a1), c2 = cross (a2 - a1, b2 - a1);
    return dcmp (c1 * c2) <= 0;
}
bool on_seg(Point p, Point a, Point b)    {         //判断点在线段上，两点式
    return dcmp (cross (a - p, b - p)) == 0 && dcmp (dot (a - p, b - p)) < 0;
}
double area_triangle(Point a, Point b, Point c) {       //三角形面积，叉积
    return fabs (cross (b - a, c - a)) / 2.0;
}
double area_poly(vector<Point> ps)  {       //多边形面积，叉积
    double ret = 0;
    for (int i=1; i<ps.size ()-1; ++i)  {
        ret += fabs (cross (ps[i] - ps[0], ps[i+1] - ps[0])) / 2;
    }
    return ret;
}
/*
    点集凸包，输入点的集合，返回凸包点的集合。
    凸包边上无点：<=    凸包边上有点：<
*/
vector<Point> convex_hull(vector<Point> ps) {
    sort (ps.begin (), ps.end ());      //x - y排序
    ps.erase (unique (ps.begin (), ps.end ()), ps.end ());  //删除重复点
    int n = ps.size (), k = 0;
    vector<Point> qs (n * 2);
    for (int i=0; i<n; ++i) {
        while (k > 1 && cross (qs[k-1] - qs[k-2], ps[i] - qs[k-2]) <= 0)  k--;
        qs[k++] = ps[i];
    }
    for (int t=k, i=n-2; i>=0; --i)  {
        while (k > t && cross (qs[k-1] - qs[k-2], ps[i] - qs[k-2]) <= 0)  k--;
        qs[k++] = ps[i];
    }
    qs.resize (k-1);
    return qs;
}
/*
	*最近点对问题，分治法，先按照x坐标排序，求解(left, mid)和(mid+1, right)范围的最小值，
	*然后类似区间合并，分离mid左右的点也求最小值
*/
double min_dist(int left, int right) {
    if (left == right) {
        return INF;
    }
    else if (right - left == 1) {
        return get_dist (point[left], point[right]);
    } else {
        int mid = left + right >> 1;
        double ret = std::min (min_dist (left, mid), min_dist (mid + 1, right));
        if (ret == 0) {
            return ret;
        }
        int endy = 0;
        for (int i=mid; i>=left&&point[mid].x-point[i].x<=ret; --i) {
            idy[endy++] = i;
        }
        for (int i=mid+1; i<=right&&point[i].x-point[mid+1].x<=ret; ++i) {
            idy[endy++] = i;
        }
        std::sort (idy, idy+endy, cmp_y);  //return point[i].y < point[j].y;
        for (int i=0; i<endy; ++i) {
            for (int j=i+1; j<endy&&point[j].y-point[i].y<ret; ++j) {
                ret = std::min (ret, get_dist (point[i], point[j]));
            }
        }
        return ret;
    }
}

struct Line {
    Point p;
    Vector v;
    double ang;
    Line () {}
    Line (const Point &p, const Vector &v) : p (p), v (v)   {
        ang = polar_angle (v);
    }
    bool operator < (const Line &r) const {
        return ang < r.ang;
    }
    Point point(double a)   {
        return p + v * a;
    }
};
/*
    半平面交，输入直线集，返回凸多边形的点集，如果不存在返回的是空集
*/
vector<Point> half_plane_inter(vector<Line> L)    {
    sort (L.begin (), L.end ());
    int first, last, n = L.size ();
    Point *p = new Point[n];
    Line *q = new Line[n];
    q[first=last=0] = L[0];
    for (int i=1; i<n; ++i) {
        while (first < last && !point_on_left (p[last-1], L[i]))    last--;
        while (first < last && !point_on_left (p[first], L[i]))  first++;
        q[++last] = L[i];
        if (dcmp (cross (q[last].v, q[last-1].v)) == 0) {
            last--;
            if (point_on_left (L[i].p, q[last]))    q[last] = L[i];
        }
        if (first < last)   p[last-1] = line_line_inter (q[last-1].p, q[last-1].v, q[last].p, q[last].v);
    }
    while (first < last && !point_on_left (p[last-1], q[first]))    last--;
    vector<Point> ps;
    if (last - first <= 1)  return ps;
    p[last] = line_line_inter (q[last].p, q[last].v, q[first].p, q[first].v);
    for (int i=first; i<=last; ++i) ps.push_back (p[i]);
    return ps;
}
bool point_on_left(Point p, Line L) {
    return cross (L.v, p - L.p) > 0;
}

struct Circle   {
    Point c;
    double r;
    Circle () {}
    Circle (Point c, double r) : c (c), r (r) {}
    Point point(double a)   {
        return Point (c.x + cos (a) * r, c.y + sin (a) * r);
    }
};
/*
    直线与圆相交求交点，返回交点个数，交点保存在P中
*/
int line_cir_inter(Line L, Circle C, double &t1, double &t2, vector<Point> &P)    {
    double a = L.v.x, b = L.p.x - C.c.x, c = L.v.y, d = L.p.y - C.c.y;
    double e = a * a + c * c, f = 2 * (a * b + c * d), g = b * b + d * d - C.r * C.r;
    double delta = f * f - 4 * e * g;
    if (dcmp (delta) < 0)   return 0;
    if (dcmp (delta) == 0)  {
        t1 = t2 = -f / (2 * e); P.push_back (L.point (t1));
        return 1;
    }
    t1 = (-f - sqrt (delta)) / (2 * e);
    t2 = (-f + sqrt (delta)) / (2 * e);
    if (t1 > t2)    swap (t1, t2);
    if (dcmp (t1) > 0 && dcmp (t1 - 1) < 0) P.push_back (L.point (t1));
    if (dcmp (t2) > 0 && dcmp (t2 - 1) < 0) P.push_back (L.point (t2));
    return (int) P.size ();
}

/*
    两圆相交求交点，返回交点个数。交点保存在P中
*/
int cir_cir_inter(Circle C1, Circle C2, vector<Point> &P)    {
    double d = length (C1.c - C2.c);
    if (dcmp (d) == 0)  {
        if (dcmp (C1.r - C2.r) == 0)    return -1;      //两圆重叠
        else    return 0;
    }
    if (dcmp (C1.r + C2.r - d) < 0) return 0;
    if (dcmp (fabs (C1.r - C2.r) - d) < 0)  return 0;
    double a = polar_angle (C2.c - C1.c);
    double da = acos ((C1.r * C1.r + d * d - C2.r * C2.r) / (2 * C1.r * d));        //C1C2到C1P1的角？
    Point p1 = C1.point (a - da), p2 = C2.point (a + da);
    P.push_back (p1);
    if (p1 == p2)   return 1;
    else    P.push_back (p2);
    return 2;
}
/*
    过点到圆的切线，返回切线条数，切线保存在V中
*/
int point_cir_tan(Point p, Circle C, Vector *V) {
    Vector u = C.c - p;
    double dis = length (u);
    if (dis < C.r)  return 0;
    else if (dcmp (dis - C.r) == 0) {
        V[0] = rotate (u, PI / 2);  return 1;
    }
    else    {
        double ang = asin (C.r / dis);
        V[0] = rotate (u, -ang);
        V[1] = rotate (u, +ang);
        return 0;
    }
}
/*
    两圆的公切线，返回公切线条数，切线短点保存在a和b中
*/
int cir_cir_tan(Circle A, Circle B, Point *a, Point *b) {
    int cnt = 0;
    if (A.r < B.r)  {
        swap (A, B);    swap (a, b);
    }
    double d = dot (A.c - B.c, A.c - B.c);
    double rsub = A.r - B.r, rsum = A.r + B.r;
    if (dcmp (d - rsub) < 0)   return 0;   //内含
    double base = polar_angle (B.c - A.c);
    if (dcmp (d) == 0 && dcmp (A.r - B.r) == 0) return -1;  //两圆重叠
    if (dcmp (d - rsub) == 0)   {       //内切，一条切线
        a[cnt] = A.point (base);    b[cnt] = B.point (base);    cnt++;
        return 1;
    }
    //有外公切线
    double ang = acos (rsub / d);
    a[cnt] = A.point (base + ang);  b[cnt] = B.point (base + ang);  cnt++;
    a[cnt] = A.point (base - ang);  b[cnt] = B.point (base - ang);  cnt++;
    if (d == rsum)  {
        a[cnt] = A.point (base);    b[cnt] = B.point (base + PI);   cnt++;
    }
    else if (dcmp (d - rsum) > 0)   {       //两条内公切线
        double ang2 = acos (rsum / d);
        a[cnt] = A.point (base + ang2); b[cnt] = B.point (base + ang2 + PI);    cnt++;
        a[cnt] = A.point (base - ang2); b[cnt] = B.point (base - ang2 + PI);    cnt++;
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
    while (dcmp (theta) <= 0)   theta += 2 * PI;
    while (theta > 2 * PI)  theta -= 2 * PI;
    theta = min (theta, 2 * PI - theta);
    return r * r * theta / 2;
}
double cal(Point a, Point b, double r)    {
    double t1, t2;
    bool ina = dcmp (length (a) - r) < 0;
    bool inb = dcmp (length (b) - r) < 0;
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
        int sgn = dcmp (cross (ps[i], ps[i+1]));
        if (sgn != 0)   {
            ret += sgn * cal (ps[i], ps[i+1], r);
        }
    }
    return ret;
}
/*************************************************/
