#ifndef RK_GEOMETRY_2D_RATIONAL_HPP
#define RK_GEOMETRY_2D_RATIONAL_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <rklib/utility/utility.hpp>
#include <vector>

namespace rklib {

namespace geo2d_rational {

template <class T>
int sgn(T x) {
    if (x < 0) return -1;
    if (x > 0) return 1;
    return 0;
}

template <class T>
struct Point {
    T x, y;

    Point(T x = 0, T y = 0) : x(x), y(y) {}

    bool operator==(const Point& rhs) const { return x == rhs.x && y == rhs.y; }
    bool operator<(const Point& rhs) const {
        if (x == rhs.x) return y < rhs.y;
        return x < rhs.x;
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
    Point& operator*=(const T& rhs) {
        this->x *= rhs;
        this->y *= rhs;
        return *this;
    }
    Point& operator/=(const T& rhs) {
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
    friend Point operator*(const Point& lhs, const T& rhs) {
        return Point(lhs) *= rhs;
    }
    friend Point operator*(const T& lhs, const Point& rhs) {
        return Point(rhs) *= lhs;
    }
    friend Point operator/(const Point& lhs, const T& rhs) {
        return Point(lhs) /= rhs;
    }
    friend T operator*(const Point& lhs, const Point& rhs) {
        return lhs.x * rhs.x + lhs.y * rhs.y;
    }
    friend T operator^(const Point& lhs, const Point& rhs) {
        return lhs.x * rhs.y - lhs.y * rhs.x;
    }

    friend std::istream& operator>>(std::istream& is, Point& p) {
        is >> p.x >> p.y;
        return is;
    }
};

template <class T>
using Vec = Point<T>;
template <class T>
using Points = std::vector<Point<T>>;
template <class T>
using Polygon = Points<T>;

template <class T>
struct Line {
    Point<T> a;
    Vec<T> v;

    Line() {}
    Line(Point<T> a, Vec<T> v) : a(a), v(v) {}
};

template <class T>
struct Segment {
    Point<T> a, b;

    Segment() {}
    Segment(Point<T> a, Point<T> b) : a(a), b(b) {}
};

template <class T>
T norm(Point<T> p) {
    return p.x * p.x + p.y * p.y;
}

template <class T>
T abs_sqared(Point<T> p) {
    return norm(p);
}

template <class T>
Point<T> projection(Line<T> l, Point<T> p) {
    auto [a, v] = l;
    T t = v * (p - a) / norm(v);
    return a + v * t;
}

template <class T>
Point<T> reflection(Line<T> l, Point<T> p) {
    return projection(l, p) * 2 - p;
}

constexpr int CCW_COUNTER_CLOCKWISE = 1;
constexpr int CCW_CLOCKWISE = -1;
constexpr int CCW_ONLINE_BACK = -2;  // C->A->B
constexpr int CCW_ONLINE_FRONT = 2;  // A->B->C
constexpr int CCW_ON_SEGMENT = 0;    // A->C->B

template <class T>
int ccw(Point<T> a, Point<T> b, Point<T> c) {
    Vec<T> v = b - a, w = c - a;
    if ((v ^ w) > T(0)) return CCW_COUNTER_CLOCKWISE;
    if ((v ^ w) < T(0)) return CCW_CLOCKWISE;
    if (v * w < T(0)) return CCW_ONLINE_BACK;
    if ((a - b) * (c - b) < T(0)) return CCW_ONLINE_FRONT;
    return CCW_ON_SEGMENT;
}

template <class T>
bool is_parallel(Vec<T> v, Vec<T> w) {
    return (v ^ w) == T(0);
}

template <class T>
bool is_orthogonal(Vec<T> v, Vec<T> w) {
    return v * w == T(0);
}

template <class T>
bool has_intersection(Line<T> l, Segment<T> s) {
    auto [p, v] = l;
    auto [a, b] = s;
    return sgn(v ^ (a - p)) * sgn(v ^ (b - p)) <= 0;
}

template <class T>
bool has_intersection(Segment<T> s1, Segment<T> s2) {
    auto [a, b] = s1;
    auto [c, d] = s2;
    return ccw(a, b, c) * ccw(a, b, d) <= 0 && ccw(c, d, a) * ccw(c, d, b) <= 0;
}

template <class T>
Point<T> intersection(Line<T> l1, Line<T> l2) {
    auto [a, v] = l1;
    auto [b, w] = l2;
    T t = ((b - a) ^ w) / (v ^ w);
    return a + v * t;
}

template <class T>
Point<T> intersection(Segment<T> s1, Segment<T> s2) {
    auto [a, b] = s1;
    auto [c, d] = s2;
    return intersection(Line(a, b - a), Line(c, d - c));
}

template <class T>
T distance_squared(Line<T> l, Point<T> p) {
    return abs_sqared(projection(l, p) - p);
}

template <class T>
T distance_squared(Segment<T> s, Point<T> p) {
    auto [a, b] = s;
    if ((b - a) * (p - a) < T(0)) return abs_sqared(p - a);
    if ((a - b) * (p - b) < T(0)) return abs_sqared(p - b);
    return distance_squared(Line(a, b - a), p);
}

template <class T>
T distance_squared(Line<T> l1, Line<T> l2) {
    auto [a, v] = l1;
    auto [b, w] = l2;
    if (is_parallel(v, w)) return distance_squared(l1, b);
    return 0;
}

template <class T>
T distance_squared(Line<T> l, Segment<T> s) {
    if (has_intersection(l, s)) return 0;
    return std::min(distance_squared(l, s.a), distance_squared(l, s.b));
}

template <class T>
T distance_squared(Segment<T> s1, Segment<T> s2) {
    if (has_intersection(s1, s2)) return 0;
    return std::min({distance_squared(s1, s2.a), distance_squared(s1, s2.b),
                     distance_squared(s2, s1.a), distance_squared(s2, s1.b)});
}

template <class T>
T area(Polygon<T>& p) {
    T res = 0;
    for (size_t i = 0; i < p.size(); i++) {
        res += p[i] ^ p[(i + 1) % p.size()] / 2;
    }
    return std::abs(res);
}

template <class T>
T area_doubled(Polygon<T>& p) {
    T res = 0;
    for (size_t i = 0; i < p.size(); i++) {
        res += p[i] ^ p[(i + 1) % p.size()];
    }
    return std::abs(res);
}

template <class T>
bool is_convex(Polygon<T>& p) {
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

template <class T>
int point_in_polygon(Point<T> a, Polygon<T>& p) {
    int n = p.size(), wn = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        if (distance_squared(Segment(p[i], p[j]), a) == 0)
            return 1;
        else if (p[i].y <= a.y && a.y < p[j].y) {
            wn += (ccw(a, p[i], p[j]) == CCW_COUNTER_CLOCKWISE);
        } else if (p[j].y <= a.y && a.y < p[i].y) {
            wn -= (ccw(a, p[i], p[j]) == CCW_CLOCKWISE);
        }
    }
    return wn == 0 ? 0 : 2;
}

template <class T>
Polygon<T> convex_hull(Points<T> p) {
    int n = p.size();
    std::sort(p.begin(), p.end());
    Polygon<T> ch(2 * n);
    int k = 0;
    for (int i = 0; i < n; i++) {
        while (k > 1 && ((ch[k - 1] - ch[k - 2]) ^ (p[i] - ch[k - 1])) <= 0)
            --k;
        ch[k++] = p[i];
    }
    for (int i = n - 2, t = k; i >= 0; --i) {
        while (k > t && ((ch[k - 1] - ch[k - 2]) ^ (p[i] - ch[k - 1])) <= 0)
            --k;
        ch[k++] = p[i];
    }
    ch.resize(k - 1);
    return ch;
}

template <class T>
T convex_cut(Polygon<T>& p, Line<T> l) {
    int n = p.size();
    auto [a, v] = l;
    Polygon<T> q;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        if ((v ^ (p[i] - a)) >= 0) q.push_back(p[i]);
        if (has_intersection(l, Segment(p[i], p[j])) &&
            !is_parallel(v, p[j] - p[i])) {
            q.push_back(intersection(l, {p[i], p[j] - p[i]}));
        }
    }
    return area(q);
}

template <class T>
T convex_cut_doubled(Polygon<T>& p, Line<T> l) {
    int n = p.size();
    auto [a, v] = l;
    Polygon<T> q;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        if ((v ^ (p[i] - a)) >= 0) q.push_back(p[i]);
        if (has_intersection(l, Segment(p[i], p[j])) &&
            !is_parallel(v, p[j] - p[i])) {
            q.push_back(intersection(l, {p[i], p[j] - p[i]}));
        }
    }
    return area_doubled(q);
}

template <class T>
Point<T> circumcenter(Point<T> a, Point<T> b, Point<T> c) {
    Point<T> m = (a + b) / 2, n = (a + c) / 2;
    Vec<T> v = {-(b - a).y, (b - a).x}, w = {-(c - a).y, (c - a).x};
    return intersection(Line(m, v), Line(n, w));
}

// template <class T>
// Point<T> incenter(Point<T> a, Point<T> b, Point<T> c) {
//     T A = abs(b - c), B = abs(c - a), C = abs(a - b);
//     return (a * A + b * B + c * C) / (A + B + C);
// }

template <class T>
Point<T> orthocenter(Point<T> a, Point<T> b, Point<T> c) {
    Vec<T> v = {-(c - b).y, (c - b).x}, w = {-(c - a).y, (c - a).x};
    return intersection(Line(a, v), Line(b, w));
}

namespace internal {

template <class T>
std::pair<T, std::pair<int, int>> closest_pair_rec(
    std::vector<std::pair<Point<T>, int>>& p, int l, int r) {
    if (r - l <= 1) return {-1, {p.size(), p.size()}};

    int m = (l + r) / 2;
    T x = p[m].first.x;
    auto d = std::min(closest_pair_rec(p, l, m), closest_pair_rec(p, m, r));
    auto cmp = [](std::pair<Point<T>, int> a, std::pair<Point<T>, int> b) {
        return a.first.y < b.first.y;
    };
    std::inplace_merge(p.begin() + l, p.begin() + m, p.begin() + r, cmp);

    std::vector<std::pair<Point<T>, int>> q;
    for (int i = l; i < r; ++i) {
        if (std::abs(p[i].first.x - x) > d.first) continue;
        for (int j = int(q.size()) - 1; j >= 0; --j) {
            T dy = p[i].first.y - q[j].first.y;
            if (dy >= d.first) break;
            chmin_non_negative(
                d, {abs(p[i].first - q[j].first), {p[i].second, q[j].second}});
        }
        q.push_back(p[i]);
    }
    return d;
}

}  // namespace internal

template <class T>
std::pair<T, std::pair<int, int>> closest_pair(Points<T>& p) {
    std::vector<std::pair<Point<T>, int>> pid(p.size());
    for (size_t i = 0; i < p.size(); i++) {
        pid[i] = {p[i], i};
    }
    std::sort(pid.begin(), pid.end());
    return internal::closest_pair_rec(pid, 0, int(p.size()));
}

template <class T>
std::pair<T, std::pair<int, int>> farthest_pair(Polygon<T>& p) {
    auto ch = convex_hull(p);
    int n = ch.size();
    if (n == 2) {
        return {abs_sqared(ch[0] - ch[1]), {0, 1}};
    }

    int i = 0, j = 0;
    for (int k = 0; k < n; ++k) {
        if (ch[k].x < ch[i].x) i = k;
        if (ch[k].x > ch[j].x) j = k;
    }
    T d = 0;
    int a = i, b = j, si = i, sj = j;
    while (i != sj || j != si) {
        if (chmax(d, abs_sqared(ch[i] - ch[j]))) a = i, b = j;
        if (((ch[(i + 1) % n] - ch[i]) ^ (ch[(j + 1) % n] - ch[j])) < 0) {
            i = (i + 1) % n;
        } else {
            j = (j + 1) % n;
        }
    }
    return {d, {a, b}};
}

template <class T>
auto arg_cmp = [](Point<T>& a, Point<T>& b) -> bool {
    if (a.y < 0) {
        if (b.y < 0)
            return (a ^ b) > 0;
        else
            return true;
    } else if (a.y == 0) {
        if (b.y < 0)
            return false;
        else if (b.y == 0)
            return a.x > b.x;
        else
            return a.x >= 0;
    } else {
        if (b.y < 0)
            return false;
        else if (b.y == 0)
            return b.x < 0;
        else
            return (a ^ b) > 0;
    }
};

template <class T>
void arg_sort(Points<T>& p) {
    std::sort(p.begin(), p.end(), arg_cmp<T>);
}

template <class T>
struct Circle {
    Point<T> c;
    T r;
};

template <class T>
int num_of_intersection(Circle<T> c, Line<T> l) {
    T d = distance_squared(l, c.c);
    if (d <= c.r) return 2;
    if (d == c.r) return 1;
    return 0;
}

template <class T>
bool has_intersection(Circle<T> c, Segment<T> s) {
    auto [a, b] = s;
    return distance_squared(s, c.c) <= c.r &&
           std::max(abs_sqared(a - c.c), abs_sqared(b - c.c)) >= c.r;
}

template <class T>
Circle<T> minimum_bounding_circle(Points<T>& p) {
    if (p.size() == 1) {
        return {p[0], 0};
    } else if (p.size() == 2) {
        Point<T> c = (p[0] + p[1]) / 2;
        return {c, abs_sqared(p[0] - c)};
    } else {
        Point<T> c;
        T r = -1;
        Points<T> ch = convex_hull(p);
        int K = ch.size();
        auto check = [&](Point<T> tc, T tr) {
            for (int i = 0; i < K; ++i) {
                if (abs_sqared(ch[i] - tc) > tr) return false;
            }
            return true;
        };
        for (int i = 0; i < K; ++i) {
            for (int j = i + 1; j < K; ++j) {
                Point<T> tc = (ch[i] + ch[j]) / 2;
                T tr = abs_sqared(ch[i] - tc);
                if (check(tc, tr) && chmin_non_negative(r, tr)) c = tc;
                for (int k = j + 1; k < K; ++k) {
                    int ccw_flag = ccw(ch[i], ch[j], ch[k]);
                    if (ccw_flag != CCW_COUNTER_CLOCKWISE &&
                        ccw_flag != CCW_CLOCKWISE)
                        continue;
                    tc = circumcenter(ch[i], ch[j], ch[k]);
                    tr = abs_sqared(ch[i] - tc);
                    if (check(tc, tr) && chmin_non_negative(r, tr)) c = tc;
                }
            }
        }
        return {c, r};
    }
}

}  // namespace geo2d_rational

}  // namespace rklib

#endif  // RK_GEOMETRY_2D_RATIONAL_HPP