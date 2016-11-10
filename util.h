#pragma once

#include "function_traits.h"

template <typename T>
struct id {
    using type = T;
};

template <typename T, T V>
struct static_value {
    using type = T;
    static constexpr T value = V;
};

template <typename T>
struct find_type;

template <typename T, typename R>
struct find_type<R T::*> {
    using type = T;
};

template <bool cond, typename T>
struct const_switch {
    using type = std::conditional_t<cond, T, T&>;
};

template <bool cond, typename T>
struct const_switch<cond, T &> {
    using type = std::conditional_t<cond, const T&, T&>;
};

template <bool cond, typename T>
using const_switch_t = typename const_switch<cond, T>::type;

