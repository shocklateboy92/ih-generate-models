#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <algorithm>

template <class C, class F>
auto transform(const C &c, F f) {
    std::vector<decltype(f(c[0]))> v(c.size());
    std::transform(c.begin(), c.end(), v.begin(), f);
    return v;
}

template <class C, class F, class I>
I accumulate(const C &c, const I &i, F f) {
    return std::accumulate(c.begin(), c.end(), i, f);
}

//template <class C, typename...A>
//auto accumulate(const C &c, A... a) {
//    return std::accumulate(c.begin, c.end(), a...);
//}

#endif // UTILS_HPP

