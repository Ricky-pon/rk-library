#ifndef RK_GEOMETRY_2D_HPP
#define RK_GEOMETRY_2D_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <rklib/utility/utility.hpp>
#include <tuple>
#include <vector>

namespace rklib {

namespace geo2d {

using real_num = long double;
real_num eps = 1e-9;
constexpr real_num PI = 3.14159265358979323846264338327950;

void set_eps(real_num new_eps) { eps = new_eps; }

int sgn(real_num x) {
    if (x < -eps) return -1;
    if (x > eps) return 1;
    return 0;
}

bool eq(real_num x, real_num y) { return sgn(x - y) == 0; }

bool ge(real_num x, real_num y) { return sgn(x - y) == 1; }

bool le(real_num x, real_num y) { return sgn(x - y) == -1; }

bool geq(real_num x, real_num y) { return sgn(x - y) >= 0; }

bool leq(real_num x, real_num y) { return sgn(x - y) <= 0; }

struct Point {
    real_num x, y;

    Point(real_num x = 0, real_num y = 0) : x(x), y(y) {}

    bool operator==(const Point& p) { return eq(x, p.x) && eq(y, p.y); }
    bool operator<(const Point& p) const {
        if (eq(x, p.x)) return le(y, p.y);
        return le(x, p.x);
    }

    Point& operator+=(const Point& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        return *this;
    }
    Point& operator-=(const Point& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        return *this;
    }
    Point& operator*=(const real_num& rhs) {
        this->x *= rhs;
        this->y *= rhs;
        return *this;
    }
    Point& operator/=(const real_num& rhs) {
        this->x /= rhs;
        this->y /= rhs;
        return *this;
    }

    friend Point operator+(const Point& lhs, const Point& rhs) {
        return Point(lhs) += rhs;
    }
    friend Point operator-(const Point& lhs, const Point& rhs) {
        return Point(lhs) -= rhs;
    }
    friend Point operator*(const Point& lhs, const real_num& rhs) {
        return Point(lhs) *= rhs;
    }
    friend Point operator*(const real_num& lhs, const Point& rhs) {
        return Point(rhs) *= lhs;
    }
    friend Point operator/(const Point& lhs, const real_num& rhs) {
        return Point(lhs) /= rhs;
    }
    friend real_num operator*(const Point& lhs, const Point& rhs) {
        return lhs.x * rhs.x + lhs.y * rhs.y;
    }
    friend real_num operator^(const Point& lhs, const Point& rhs) {
        return lhs.x * rhs.y - lhs.y * rhs.x;
    }

    friend std::istream& operator>>(std::istream& is, Point& p) {
        is >> p.x >> p.y;
        return is;
    }
};

using Vec = Point;
using Points = std::vector<Point>;
using Polygon = Points;

struct Line {
    Point a;
    Vec v;

    Line() {}
    Line(Point a, Vec v) : a(a), v(v) {}
};

struct Segment {
    Point a, b;

    Segment() {}
    Segment(Point a, Point b) : a(a), b(b) {}
};

real_num norm(Point p) { return p.x * p.x + p.y * p.y; }

real_num abs(Point p) { return std::sqrt(norm(p)); }

real_num arg(Point p) { return std::atan2(p.y, p.x); }

real_num angle(Point a, Point b) { return arg({a * b, a ^ b}); }

Point rotation(Point p, real_num t) {
    return {p.x * std::cos(t) - p.y * std::sin(t),
            p.x * std::sin(t) + p.y * std::cos(t)};
}

Point projection(Line l, Point p) {
    auto [a, v] = l;
    real_num t = v * (p - a) / norm(v);
    return a + v * t;
}

Point reflection(Line l, Point p) { return projection(l, p) * 2 - p; }

constexpr int CCW_COUNTER_CLOCKWISE = 1;
constexpr int CCW_CLOCKWISE = -1;
constexpr int CCW_ONLINE_BACK = -2;  // C->A->B
constexpr int CCW_ONLINE_FRONT = 2;  // A->B->C
constexpr int CCW_ON_SEGMENT = 0;    // A->C->B

int ccw(Point a, Point b, Point c) {
    Vec v = b - a, w = c - a;
    if (ge(v ^ w, 0)) return CCW_COUNTER_CLOCKWISE;
    if (le(v ^ w, 0)) return CCW_CLOCKWISE;
    if (le(v * w, 0)) return CCW_ONLINE_BACK;
    if (le((a - b) * (c - b), 0)) return CCW_ONLINE_FRONT;
    return CCW_ON_SEGMENT;
}

bool is_parallel(Vec v, Vec w) { return eq(v ^ w, 0); }

bool is_orthogonal(Vec v, Vec w) { return eq(v * w, 0); }

bool has_intersection(Line l, Segment s) {
    auto [p, v] = l;
    auto [a, b] = s;
    return sgn(v ^ (a - p)) * sgn(v ^ (b - p)) <= 0;
}

bool has_intersection(Segment s1, Segment s2) {
    auto [a, b] = s1;
    auto [c, d] = s2;
    return ccw(a, b, c) * ccw(a, b, d) <= 0 && ccw(c, d, a) * ccw(c, d, b) <= 0;
}

Point intersection(Line l1, Line l2) {
    auto [a, v] = l1;
    auto [b, w] = l2;
    real_num t = ((b - a) ^ w) / (v ^ w);
    return a + v * t;
}

Point intersection(Segment s1, Segment s2) {
    auto [a, b] = s1;
    auto [c, d] = s2;
    return intersection(Line(a, b - a), Line(c, d - c));
}

real_num distance(Line l, Point p) {
    auto [a, v] = l;
    return std::abs(v ^ (p - a) / abs(v));
}

real_num distance(Segment s, Point p) {
    auto [a, b] = s;
    if (le((b - a) * (p - a), 0)) return abs(p - a);
    if (le((a - b) * (p - b), 0)) return abs(p - b);
    return distance(Line(a, b - a), p);
}

real_num distance(Line l1, Line l2) {
    auto [a, v] = l1;
    auto [b, w] = l2;
    if (is_parallel(v, w)) return distance(l1, b);
    return 0;
}

real_num distance(Line l, Segment s) {
    if (has_intersection(l, s)) return 0;
    return std::min(distance(l, s.a), distance(l, s.b));
}

real_num distance(Segment s1, Segment s2) {
    if (has_intersection(s1, s2)) return 0;
    return std::min({distance(s1, s2.a), distance(s1, s2.b), distance(s2, s1.a),
                     distance(s2, s1.b)});
}

real_num area(Polygon& p) {
    real_num res = 0;
    for (size_t i = 0; i < p.size(); i++) {
        res += p[i] ^ p[(i + 1) % p.size()] / 2;
    }
    return std::abs(res);
}

bool is_convex(Polygon& p) {
    int n = p.size();
    bool flag1 = false, flag2 = false;
    for (int i = 0; i < n; i++) {
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

int point_in_polygon(Point a, Polygon& p) {
    int n = p.size(), wn = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        if (eq(distance(Segment(p[i], p[j]), a), 0)) {
            return 1;
        } else if (p[i].y <= a.y && a.y < p[j].y) {
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
    for (int i = 0; i < n; i++) {
        while (k > 1 && le((ch[k - 1] - ch[k - 2]) ^ (p[i] - ch[k - 1]), 0)) {
            --k;
        }
        ch[k++] = p[i];
    }
    for (int i = n - 2, t = k; i >= 0; --i) {
        while (k > t && le((ch[k - 1] - ch[k - 2]) ^ (p[i] - ch[k - 1]), 0)) {
            --k;
        }
        ch[k++] = p[i];
    }
    ch.resize(k - 1);
    return ch;
}

real_num convex_cut(Polygon& p, Line l) {
    int n = p.size();
    auto [a, v] = l;
    Polygon q;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        if (geq(v ^ (p[i] - a), 0)) q.push_back(p[i]);
        if (has_intersection(l, Segment(p[i], p[j])) &&
            !is_parallel(v, p[j] - p[i])) {
            q.push_back(intersection(l, {p[i], p[j] - p[i]}));
        }
    }
    return area(q);
}

Point circumcenter(Point a, Point b, Point c) {
    Point m = (a + b) / 2, n = (a + c) / 2;
    Vec v = {-(b - a).y, (b - a).x}, w = {-(c - a).y, (c - a).x};
    return intersection(Line(m, v), Line(n, w));
}

Point incenter(Point a, Point b, Point c) {
    real_num A = abs(b - c), B = abs(c - a), C = abs(a - b);
    return (a * A + b * B + c * C) / (A + B + C);
}

Point orthocenter(Point a, Point b, Point c) {
    Vec v = {-(c - b).y, (c - b).x}, w = {-(c - a).y, (c - a).x};
    return intersection(Line(a, v), Line(b, w));
}

namespace internal {

std::pair<real_num, std::pair<int, int>> closest_pair_rec(
    std::vector<std::pair<Point, int>>& p, int l, int r) {
    if (r - l <= 1)
        return {std::numeric_limits<real_num>::max(), {p.size(), p.size()}};

    int m = (l + r) / 2;
    real_num x = p[m].first.x;
    auto d = std::min(closest_pair_rec(p, l, m), closest_pair_rec(p, m, r));
    auto cmp = [](std::pair<Point, int> a, std::pair<Point, int> b) {
        return a.first.y < b.first.y;
    };
    std::inplace_merge(p.begin() + l, p.begin() + m, p.begin() + r, cmp);

    std::vector<std::pair<Point, int>> q;
    for (int i = l; i < r; ++i) {
        if (ge(std::abs(p[i].first.x - x), d.first)) continue;
        for (int j = int(q.size()) - 1; j >= 0; --j) {
            real_num dy = p[i].first.y - q[j].first.y;
            if (geq(dy, d.first)) break;
            chmin(d,
                  {abs(p[i].first - q[j].first), {p[i].second, q[j].second}});
        }
        q.push_back(p[i]);
    }
    return d;
}

}  // namespace internal

std::pair<real_num, std::pair<int, int>> closest_pair(Points& p) {
    std::vector<std::pair<Point, int>> pid(p.size());
    for (size_t i = 0; i < p.size(); i++) {
        pid[i] = {p[i], i};
    }
    std::sort(pid.begin(), pid.end());
    return internal::closest_pair_rec(pid, 0, int(p.size()));
}

std::pair<real_num, std::pair<int, int>> farthest_pair(Polygon& p) {
    int n = p.size();
    if (n == 2) {
        return {abs(p[0] - p[1]), {0, 1}};
    }
    int i = 0, j = 0;
    for (int k = 0; k < n; ++k) {
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

void arg_sort(Points& p) {
    auto cmp = [&](Point& a, Point& b) {
        if (le(a.y, 0)) {
            if (le(b.y, 0)) {
                return ge(a ^ b, 0);
            } else {
                return true;
            }
        } else if (eq(a.y, 0)) {
            if (le(b.y, 0)) {
                return false;
            } else if (eq(b.y, 0)) {
                return ge(a.x, b.x);
            } else {
                return geq(a.x, 0);
            }
        } else {
            if (le(b.y, 0)) {
                return false;
            } else if (eq(b.y, 0)) {
                return le(b.x, 0);
            } else {
                return ge(a ^ b, 0);
            }
        }
    };
    std::sort(p.begin(), p.end(), cmp);
}

struct Circle {
    Point c;
    real_num r;
};

int num_of_common_tangent(Circle c1, Circle c2) {
    if (c1.r < c2.r) std::swap(c1, c2);
    real_num d = abs(c1.c - c2.c), r = c1.r + c2.r;
    if (ge(d, r)) return 4;
    if (eq(d, r)) return 3;
    if (eq(d + c2.r, c1.r)) return 1;
    if (le(d + c2.r, c1.r)) return 0;
    return 2;
}

int num_of_intersection(Circle c, Line l) {
    real_num d = distance(l, c.c);
    if (leq(d, c.r)) return 2;
    if (eq(d, c.r)) return 1;
    return 0;
}

bool has_intersection(Circle c, Segment s) {
    auto [a, b] = s;
    return leq(distance(s, c.c), c.r) &&
           geq(std::max(abs(a - c.c), abs(b - c.c)), c.r);
}

Points intersection(Circle c, Line l) {
    Points res;
    if (num_of_intersection(c, l) == 0) return res;
    Point p = projection(l, c.c);
    real_num t = std::sqrt(
        std::max(real_num(0), (c.r * c.r - norm(p - c.c)) / norm(l.v)));
    res.push_back(p + l.v * t);
    if (!eq(t, 0)) res.push_back(p - l.v * t);
    return res;
}

Points intersection(Circle c, Segment s) {
    auto [a, b] = s;
    Points ps = intersection(c, Line(b, a - b));
    Points res;
    for (auto p : ps) {
        if (ccw(a, b, p) == CCW_ON_SEGMENT) res.push_back(p);
    }
    return res;
}

Points intersection(Circle c1, Circle c2) {
    real_num r1 = c1.r, r2 = c2.r;
    Points res;
    Vec v = c2.c - c1.c, w = {-v.y, v.x};
    real_num d = abs(v);
    real_num x = (d * d + r1 * r1 - r2 * r2) / (2 * d);
    real_num y = std::sqrt(std::max(r1 * r1 - x * x, real_num(0)));
    res.push_back(c1.c + v * x / d + w * y / d);
    if (num_of_common_tangent(c1, c2) != 2) return res;
    res.push_back(c1.c + v * x / d - w * y / d);
    return res;
}

namespace internal {

real_num area_of_intersection(Circle c, Segment s) {
    auto [a, b] = s;
    Vec va = a - c.c, vb = b - c.c;
    if (eq(va ^ vb, 0))
        return 0;
    else if (leq(abs(va), c.r) && leq(abs(vb), c.r))
        return (va ^ vb) / 2;
    else if (geq(distance(s, c.c), c.r))
        return c.r * c.r * angle(va, vb) / 2;
    else {
        auto ps = intersection(c, s);
        if (ps.size() == 1)
            return area_of_intersection(c, {a, ps[0]}) +
                   area_of_intersection(c, {ps[0], b});
        else
            return area_of_intersection(c, {a, ps[0]}) +
                   area_of_intersection(c, {ps[0], ps[1]}) +
                   area_of_intersection(c, {ps[1], b});
    }
}

}  // namespace internal

real_num area_of_intersection(Circle c, Polygon& p) {
    int n = p.size();
    real_num res = 0;
    for (int i = 0; i < n; ++i) {
        res += internal::area_of_intersection(c, {p[i], p[(i + 1) % n]});
    }
    return res;
}

real_num area_of_intersection(Circle c1, Circle c2) {
    int num = num_of_common_tangent(c1, c2);
    if (num >= 3) return 0;
    real_num r1 = c1.r, r2 = c2.r;
    if (num <= 1) {
        real_num r = std::min(r1, r2);
        return PI * r * r;
    }
    real_num d = abs(c1.c - c2.c);
    real_num res = 0;
    for (int i = 0; i < 2; ++i) {
        real_num x = (d * d + r1 * r1 - r2 * r2) / (2 * d);
        real_num t = std::acos(x / r1) * 2;
        res += (t - std::sin(t)) * r1 * r1 / 2;
        std::swap(r1, r2);
    }
    return res;
}

Points tangent(Circle c, Point p) {
    Points res;
    real_num d = abs(p - c.c);
    real_num t = std::acos(c.r / d);
    res.push_back(c.c + rotation(p - c.c, t) * c.r / d);
    res.push_back(c.c + rotation(p - c.c, -t) * c.r / d);
    return res;
}

std::tuple<std::vector<Line>, Points, Points> common_tangent(Circle c1,
                                                             Circle c2) {
    std::vector<Line> ls;
    Points ps, qs;
    int num = num_of_common_tangent(c1, c2);
    if (num >= 2) {
        real_num d = abs(c2.c - c1.c);
        real_num t = std::acos(std::abs(c1.r - c2.r) / d);
        if (le(c1.r, c2.r)) t = PI - t;

        ps.push_back(c1.c + rotation(c2.c - c1.c, t) * c1.r / d);
        qs.push_back(c2.c + rotation(c2.c - c1.c, t) * c2.r / d);
        ls.push_back({ps.back(), qs.back() - ps.back()});

        ps.push_back(c1.c + rotation(c2.c - c1.c, -t) * c1.r / d);
        qs.push_back(c2.c + rotation(c2.c - c1.c, -t) * c2.r / d);
        ls.push_back({ps.back(), qs.back() - ps.back()});
    }
    if (num == 4) {
        real_num d = abs(c2.c - c1.c);
        real_num L = d * c1.r / (c1.r + c2.r);
        real_num t = std::acos(c1.r / L);

        ps.push_back(c1.c + rotation(c2.c - c1.c, t) * c1.r / d);
        qs.push_back(c2.c + rotation(c1.c - c2.c, t) * c2.r / d);
        ls.push_back({ps.back(), qs.back() - ps.back()});

        ps.push_back(c1.c + rotation(c2.c - c1.c, -t) * c1.r / d);
        qs.push_back(c2.c + rotation(c1.c - c2.c, -t) * c2.r / d);
        ls.push_back({ps.back(), qs.back() - ps.back()});
    }
    if (num == 3 || num == 1) {
        Points tg = intersection(c1, c2);
        ps.push_back(tg[0]);
        qs.push_back(tg[0]);
        ls.push_back({tg[0], rotation(c2.c - c1.c, PI / 2)});
    }
    return {ls, ps, qs};
}

Circle minimum_bounding_circle(Points& p) {
    if (p.size() == 1) {
        return {p[0], 0};
    } else if (p.size() == 2) {
        Point c = (p[0] + p[1]) / 2;
        return {c, abs(p[0] - c)};
    } else {
        Point c;
        real_num r = std::numeric_limits<real_num>::max();
        Points ch = convex_hull(p);
        int K = ch.size();
        auto check = [&](Point tc, real_num tr) {
            for (int i = 0; i < K; ++i) {
                if (ge(abs(ch[i] - tc), tr)) return false;
            }
            return true;
        };
        for (int i = 0; i < K; ++i) {
            for (int j = i + 1; j < K; ++j) {
                Point tc = (ch[i] + ch[j]) / 2;
                real_num tr = abs(ch[i] - tc);
                if (check(tc, tr) && chmin(r, tr)) c = tc;
                for (int k = j + 1; k < K; ++k) {
                    int ccw_flag = ccw(ch[i], ch[j], ch[k]);
                    if (ccw_flag != CCW_COUNTER_CLOCKWISE &&
                        ccw_flag != CCW_CLOCKWISE)
                        continue;
                    tc = circumcenter(ch[i], ch[j], ch[k]);
                    tr = abs(ch[i] - tc);
                    if (check(tc, tr) && chmin(r, tr)) c = tc;
                }
            }
        }
        return {c, r};
    }
}

}  // namespace geo2d

}  // namespace rklib

#endif  // RK_GEOMETRY_2D_HPP