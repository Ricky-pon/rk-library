namespace geo {
using real_num = double;
constexpr real_num eps = 1e-9;
constexpr real_num PI = 3.14159265358979323846264338327950;

inline int sgn(real_num x) {
    if (x < -eps) return -1;
    if (x > eps) return 1;
    return 0;
}

inline bool eq(real_num x, real_num y) { return sgn(x - y) == 0; }

inline bool ge(real_num x, real_num y) { return sgn(x - y) == 1; }

inline bool le(real_num x, real_num y) { return sgn(x - y) == -1; }

inline bool geq(real_num x, real_num y) { return sgn(x - y) >= 0; }

inline bool leq(real_num x, real_num y) { return sgn(x - y) <= 0; }

struct Point {
    real_num x, y;
    Point(real_num x = 0, real_num y = 0) : x(x), y(y) {}

    Point operator+(const Point &p) { return {x + p.x, y + p.y}; }

    Point operator-(const Point &p) { return {x - p.x, y - p.y}; }

    Point operator*(const real_num k) { return {k * x, k * y}; }

    Point operator/(const real_num k) { return {x / k, y / k}; }

    real_num operator*(const Point &p) { return x * p.x + y * p.y; }

    real_num operator^(const Point &p) { return x * p.y - y * p.x; }

    bool operator==(const Point &p) { return eq(x, p.x) && eq(y, p.y); }

    bool operator<(const Point &p) const {
        if (eq(x, p.x)) return le(y, p.y);
        return le(x, p.x);
    }
};

using Vec = Point;
using Points = vector<Point>;
using Polygon = vector<Point>;

real_num norm(Point p) { return p.x * p.x + p.y * p.y; }

real_num abs(Point p) { return sqrt(norm(p)); }

real_num arg(Point p) { return atan2(p.y, p.x); }

real_num angle(Point a, Point b) { return arg({a * b, a ^ b}); }

Point rot(Point p, real_num t) {
    return {p.x * cos(t) - p.y * sin(t), p.x * sin(t) + p.y * cos(t)};
}

Point proj(Point a, Vec v, Point p) {
    real_num t = v * (p - a) / norm(v);
    return a + v * t;
}

Point refl(Point a, Vec v, Point p) { return proj(a, v, p) * 2 - p; }

constexpr int CCW_COUNTER_CLOCKWISE = 1;
constexpr int CCW_CLOCKWISE = -1;
constexpr int CCW_ONLINE_BACK = -2;  // C->A->B
constexpr int CCW_ONLINE_FRONT = 2;  // A->B->C
constexpr int CCW_ON_SEGMENT = 0;    // A->C->B

inline int ccw(Point a, Point b, Point c) {
    Vec v = b - a, w = c - a;
    if (ge(v ^ w, 0)) return CCW_COUNTER_CLOCKWISE;
    if (le(v ^ w, 0)) return CCW_CLOCKWISE;
    if (le(v * w, 0)) return CCW_ONLINE_BACK;
    if (le((a - b) * (c - b), 0)) return CCW_ONLINE_FRONT;
    return CCW_ON_SEGMENT;
}

bool is_parallel(Vec v, Vec w) { return eq(v ^ w, 0); }

bool is_orthogonal(Vec v, Vec w) { return eq(v * w, 0); }

bool has_intersection_ls(Point p, Vec v, Point a, Point b) {
    return sgn(v ^ (a - p)) * sgn(v ^ (b - p)) <= 0;
}

bool has_intersection_ss(Point a, Point b, Point c, Point d) {
    return ccw(a, b, c) * ccw(a, b, d) <= 0 && ccw(c, d, a) * ccw(c, d, b) <= 0;
}

Point intersection_ll(Point a, Vec v, Point b, Vec w) {
    real_num t = ((b - a) ^ w) / (v ^ w);
    return a + v * t;
}

real_num distance_lp(Point a, Vec v, Point p) {
    return abs(v ^ (p - a) / abs(v));
}

real_num distance_sp(Point a, Point b, Point p) {
    if (le((b - a) * (p - a), 0)) return abs(p - a);
    if (le((a - b) * (p - b), 0)) return abs(p - b);
    return distance_lp(a, b - a, p);
}

real_num distance_ll(Point a, Vec v, Point b, Vec w) {
    if (is_parallel(v, w)) return distance_lp(a, v, b);
    return 0;
}

real_num distance_ls(Point p, Vec v, Point a, Point b) {
    if (has_intersection_ls(p, v, a, b)) return 0;
    return min(distance_lp(p, v, a), distance_lp(p, v, b));
}

real_num distance_ss(Point a, Point b, Point c, Point d) {
    if (has_intersection_ss(a, b, c, d)) return 0;
    return min({distance_sp(a, b, c), distance_sp(a, b, d),
                distance_sp(c, d, a), distance_sp(c, d, b)});
}

real_num area(Polygon &p) {
    real_num ret = 0;
    rep(i, p.size()) ret += p[i] ^ p[(i + 1) % p.size()] / 2;
    return abs(ret);
}

bool is_convex(Polygon &p) {
    int n = p.size();
    bool flag1 = false, flag2 = false;
    rep(i, n) {
        int tmp = ccw(p[(i + n - 1) % n], p[i], p[(i + 1) % n]);
        if (tmp == CCW_COUNTER_CLOCKWISE) {
            if (flag2) return false;
            flag1 = true;
        } else if (tmp == CCW_CLOCKWISE) {
            if (flag1) return false;
            flag2 = true;
        }
    }
    return true;
}

int point_in_polygon(Point a, Polygon &p) {
    int n = p.size(), wn = 0;
    rep(i, n) {
        int j = (i + 1) % n;
        if (distance_sp(p[i], p[j], a) == 0)
            return 1;
        else if (p[i].y <= a.y && a.y < p[j].y) {
            wn += (ccw(a, p[i], p[j]) == CCW_COUNTER_CLOCKWISE);
        } else if (p[j].y <= a.y && a.y < p[i].y) {
            wn -= (ccw(a, p[i], p[j]) == CCW_CLOCKWISE);
        }
    }
    return wn == 0 ? 0 : 2;
}

Polygon convex_hull(Points p) {
    int n = p.size();
    sort(p.begin(), p.end());
    Polygon ch(2 * n);
    int k = 0;
    rep(i, n) {
        while (k > 1 && le((ch[k - 1] - ch[k - 2]) ^ (p[i] - ch[k - 1]), 0))
            --k;
        ch[k++] = p[i];
    }
    for (int i = n - 2, t = k; i >= 0; --i) {
        while (k > t && le((ch[k - 1] - ch[k - 2]) ^ (p[i] - ch[k - 1]), 0))
            --k;
        ch[k++] = p[i];
    }
    ch.resize(k - 1);
    return ch;
}

pair<real_num, pii> farthest_pair(Polygon &p) {
    int n = p.size();
    if (n == 2) {
        return {abs(p[0] - p[1]), {0, 1}};
    }
    int i = 0, j = 0;
    rep(k, n) {
        if (le(p[k].x, p[i].x)) i = k;
        if (ge(p[k].x, p[j].x)) j = k;
    }
    real_num d = 0;
    int a = i, b = j, si = i, sj = j;
    while (i != sj || j != si) {
        if (chmax(d, abs(p[i] - p[j]))) a = i, b = j;
        if (le((p[(i + 1) % n] - p[i]) ^ (p[(j + 1) % n] - p[j]), 0)) {
            i = (i + 1) % n;
        } else
            j = (j + 1) % n;
    }
    return {d, {a, b}};
}

real_num convex_cut(Polygon &p, Point a, Vec v) {
    int n = p.size();
    Polygon q;
    rep(i, n) {
        int j = (i + 1) % n;
        if (geq(v ^ (p[i] - a), 0)) q.push_back(p[i]);
        if (has_intersection_ls(a, v, p[i], p[j]) &&
            !is_parallel(v, p[j] - p[i])) {
            q.push_back(intersection_ll(a, v, p[i], p[j] - p[i]));
        }
    }
    return area(q);
}

pair<real_num, pii> closest_pair_rec(vector<pair<Point, int>> &p, int l,
                                     int r) {
    if (r - l <= 1) return {INF, {p.size(), p.size()}};

    int m = (l + r) / 2;
    real_num x = p[m].fi.x;
    auto d = min(closest_pair_rec(p, l, m), closest_pair_rec(p, m, r));
    auto cmp = [](pair<Point, int> a, pair<Point, int> b) {
        return a.fi.y < b.fi.y;
    };
    inplace_merge(p.begin() + l, p.begin() + m, p.begin() + r, cmp);

    vector<pair<Point, int>> q;
    For(i, l, r) {
        if (ge(abs(p[i].fi.x - x), d.fi)) continue;
        rrep(j, q.size()) {
            real_num dy = p[i].fi.y - q[j].fi.y;
            if (geq(dy, d.fi)) break;
            chmin(d, {abs(p[i].fi - q[j].fi), {p[i].se, q[j].se}});
        }
        q.push_back(p[i]);
    }
    return d;
}

pair<real_num, pii> closest_pair(Points &p) {
    vector<pair<Point, int>> pid(p.size());
    rep(i, p.size()) pid[i] = {p[i], i};
    sort(pid.begin(), pid.end());
    return closest_pair_rec(pid, 0, p.size());
}

int has_intersection_cc(Point c1, real_num r1, Point c2, real_num r2) {
    if (r1 < r2) {
        swap(c1, c2);
        swap(r1, r2);
    }
    real_num d = abs(c1 - c2), r = r1 + r2;
    if (ge(d, r)) return 4;
    if (eq(d, r)) return 3;
    if (eq(d + r2, r1)) return 1;
    if (le(d + r2, r1)) return 0;
    return 2;
}

bool has_intersection_cl(Point c, real_num r, Point a, Vec v) {
    return leq(distance_lp(a, v, c), r);
}

bool has_intersection_cs(Point c, real_num r, Point a, Point b) {
    return leq(distance_sp(a, b, c), r) && geq(max(abs(a - c), abs(b - c)), r);
}

Points intersection_cl(Point c, real_num r, Point a, Vec v) {
    Points ps;
    if (!has_intersection_cl(c, r, a, v)) return ps;
    Point p = proj(a, v, c);
    real_num t = sqrt(max((real_num)0.0, (r * r - norm(p - c)) / norm(v)));
    ps.push_back(p + v * t);
    if (!eq(t, 0)) ps.push_back(p - v * t);
    return ps;
}

Points intersection_cs(Point c, real_num r, Point a, Point b) {
    Points ps = intersection_cl(c, r, b, a - b);
    Points qs;
    for (auto p : ps) {
        if (ccw(a, b, p) == CCW_ON_SEGMENT) qs.push_back(p);
    }
    return qs;
}

Points intersection_cc(Point c1, real_num r1, Point c2, real_num r2) {
    Points ps;
    Vec v = c2 - c1, w = {v.y * -1, v.x};
    real_num d = abs(v);
    real_num x = (d * d + r1 * r1 - r2 * r2) / (2 * d);
    real_num y = sqrt(max(r1 * r1 - x * x, (real_num)0.0));
    ps.push_back(c1 + v * x / d + w * y / d);
    if (has_intersection_cc(c1, r1, c2, r2) != 2) return ps;
    ps.push_back(c1 + v * x / d - w * y / d);
    return ps;
}

real_num common_area_ct(Point c, real_num r, Point a, Point b) {
    Vec va = a - c, vb = b - c;
    if (eq(va ^ vb, 0))
        return 0;
    else if (leq(abs(va), r) && leq(abs(vb), r))
        return (va ^ vb) / 2;
    else if (geq(distance_sp(a, b, c), r))
        return r * r * angle(va, vb) / 2;
    else {
        auto ps = intersection_cs(c, r, a, b);
        if (ps.size() == 1)
            return common_area_ct(c, r, a, ps[0]) +
                   common_area_ct(c, r, ps[0], b);
        else
            return common_area_ct(c, r, a, ps[0]) +
                   common_area_ct(c, r, ps[0], ps[1]) +
                   common_area_ct(c, r, ps[1], b);
    }
}

real_num common_area_cp(Point c, real_num r, Polygon &p) {
    int n = p.size();
    real_num ret = 0;
    rep(i, n) { ret += common_area_ct(c, r, p[i], p[(i + 1) % n]); }
    return ret;
}

real_num commn_area_cc(Point c1, real_num r1, Point c2, real_num r2) {
    int flag = has_intersection_cc(c1, r1, c2, r2);
    if (flag >= 3) return 0;
    if (flag <= 1) {
        real_num r = min(r1, r2);
        return PI * r * r;
    }
    real_num d = abs(c1 - c2);
    real_num ret = 0;
    rep(i, 2) {
        real_num x = (d * d + r1 * r1 - r2 * r2) / (2 * d);
        real_num t = acos(x / r1) * 2;
        ret += (t - sin(t)) * r1 * r1 / 2;
        swap(c1, c2);
        swap(r1, r2);
    }
    return ret;
}

Points tangent(Point c, real_num r, Point p) {
    Points ps;
    real_num d = abs(p - c);
    real_num t = acos(r / d);
    ps.push_back(c + rot(p - c, t) * r / d);
    ps.push_back(c + rot(p - c, -t) * r / d);
    return ps;
}

Points common_tangent(Point c1, real_num r1, Point c2, real_num r2) {
    Points ps;
    int flag = has_intersection_cc(c1, r1, c2, r2);
    if (flag >= 2) {
        real_num d = abs(c2 - c1);
        real_num t = acos(abs(r1 - r2) / d);
        if (le(r1, r2)) t = PI - t;
        ps.push_back(c1 + rot(c2 - c1, t) * r1 / d);
        ps.push_back(c1 + rot(c2 - c1, -t) * r1 / d);
    }
    if (flag == 4) {
        real_num d = abs(c2 - c1);
        real_num L = d * r1 / (r1 + r2);
        real_num t = acos(r1 / L);
        ps.push_back(c1 + rot(c2 - c1, t) * r1 / d);
        ps.push_back(c1 + rot(c2 - c1, -t) * r1 / d);
    }
    if (flag == 3 || flag == 1) {
        Polygon tg = intersection_cc(c1, r1, c2, r2);
        ps.push_back(tg[0]);
    }
    return ps;
}

Point get_o(Point a, Point b, Point c) {
    Point M = (a + b) / 2, N = (a + c) / 2;
    Vec v = {-(b - a).y, (b - a).x}, w = {-(c - a).y, (c - a).x};
    return intersection_ll(M, v, N, w);
}

Point get_i(Point a, Point b, Point c) {
    real_num A = abs(b - c), B = abs(c - a), C = abs(a - b);
    return (a * A + b * B + c * C) / (A + B + C);
}

Point get_h(Point a, Point b, Point c) {
    Vec v = {-(c - b).y, (c - b).x}, w = {-(c - a).y, (c - a).x};
    return intersection_ll(a, v, b, w);
}

pair<Point, real_num> minimum_bounding_circle(Points &p) {
    Point C;
    real_num r;
    if (p.size() == 1)
        C = p[0], r = 0;
    else if (p.size() == 2)
        C = (p[0] + p[1]) / 2, r = abs(p[0] - C);
    else {
        r = INF;
        Points ch = convex_hull(p);
        int K = ch.size();
        auto check = [&](Point tc, real_num tr) {
            rep(i, K) {
                if (ge(abs(ch[i] - tc), tr)) return false;
            }
            return true;
        };
        rep(i, K) For(j, i + 1, K) {
            Point tc = (ch[i] + ch[j]) / 2;
            real_num tr = abs(ch[i] - tc);
            if (check(tc, tr) && chmin(r, tr)) C = tc;
            For(k, j + 1, K) {
                int ccw_flag = ccw(ch[i], ch[j], ch[k]);
                if (ccw_flag != CCW_COUNTER_CLOCKWISE &&
                    ccw_flag != CCW_CLOCKWISE)
                    continue;
                tc = get_o(ch[i], ch[j], ch[k]);
                tr = abs(ch[i] - tc);
                if (check(tc, tr) && chmin(r, tr)) C = tc;
            }
        }
    }
    return {C, r};
}

void arg_sort(Points &p) {
    auto cmp = [&](Point &a, Point &b) {
        if (le(a.y, 0)) {
            if (le(b.y, 0))
                return geq(a ^ b, 0);
            else
                return true;
        } else if (eq(a.y, 0)) {
            if (le(b.y, 0))
                return false;
            else if (eq(b.y, 0))
                return geq(a.x, b.x);
            else
                return geq(a.x, 0);
        } else {
            if (le(b.y, 0))
                return false;
            else if (eq(b.y, 0))
                return le(b.x, 0);
            else
                return geq(a ^ b, 0);
        }
    };
    sort(p.begin(), p.end(), cmp);
}
}  // namespace geo