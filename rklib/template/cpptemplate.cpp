#include <bits/stdc++.h>

#define REP_OVERLOAD(arg1, arg2, arg3, arg4, NAME, ...) NAME
#define REP3(i, l, r, s)                                                       \
    for (int i = int(l), rep3_r = int(r), rep3_s = int(s); i < rep3_r;         \
         i += rep3_s)
#define REP2(i, l, r) REP3(i, l, r, 1)
#define REP1(i, n) REP2(i, 0, n)
#define rep(...) REP_OVERLOAD(__VA_ARGS__, REP3, REP2, REP1, )(__VA_ARGS__)
#define repin(i, l, r) for (int i = int(l), repin_r = int(r); i <= repin_r; ++i)

#define RREP_OVERLOAD(arg1, arg2, arg3, arg4, NAME, ...) NAME
#define RREP3(i, l, r, s)                                                      \
    for (int i = int(r) - 1, rrep3_l = int(l), rrep3_s = int(s); i >= rrep3_l; \
         i -= rrep3_s)
#define RREP2(i, l, r) RREP3(i, l, r, 1)
#define RREP1(i, n) RREP2(i, 0, n)
#define rrep(...) RREP_OVERLOAD(__VA_ARGS__, RREP3, RREP2, RREP1, )(__VA_ARGS__)
#define rrepin(i, l, r)                                                        \
    for (int i = int(r), rrepin_l = int(l); i >= rrepin_l; --i)

#define all(v) v.begin(), v.end()

#include "rklib/utility/utility.hpp"

using namespace std;
using namespace rklib;

using lint = long long;
using ulint = unsigned long long;
using pii = pair<int, int>;
using pll = pair<lint, lint>;
template<class T>
using vec = vector<T>;
template<class T>
using vvec = vec<vec<T>>;
template<class T>
using vvvec = vec<vvec<T>>;