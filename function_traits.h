#pragma once
#include <type_traits>
#include <utility>

template <typename T>
struct remove_const_reference {
    using type = T;
};

template <typename T>
struct remove_const_reference<const T &> {
    using type = T;
};

template <typename T>
struct remove_const_reference<T &&> {
    using type = T;
};

template <typename T>
struct remove_const_reference<const T &&> {
    using type = T;
};

template <typename F>
struct function_traits
    : public function_traits<std::remove_cv_t<decltype(&F::operator())>> {
};

// function pointer
template <typename R, typename C, typename... Args>
struct function_traits<R (C::*)(Args...)> : public function_traits<R(Args...)> {
};

template <typename R, typename C, typename... Args>
struct function_traits<R (C::*)(Args...) const> : public function_traits<R(Args...)> {
};

template <typename R, typename... Args>
struct function_traits<R (*)(Args...)> : public function_traits<R(Args...)> {
};

template <typename R, typename C>
struct function_traits<R (C::*)> : public function_traits<R()> {
};

template <typename R, typename... Args>
struct function_traits<R(Args...)> {
    using prototype =
        typename remove_const_reference<R>::type(typename remove_const_reference<Args>::type...);
    using return_type = typename remove_const_reference<R>::type;

    static constexpr std::size_t arity = sizeof...(Args);

    template <std::size_t N>
    struct argument {
        static_assert(N < arity, "error: invalid parameter index.");
        using type = typename std::tuple_element<N, std::tuple<Args...>>::type;
        using pure_type = typename remove_const_reference<type>::type;
    };
};

template <typename F>
using function_prototype = typename function_traits<F>::prototype;

template <typename F, typename R>
struct same_prototype {
    static constexpr bool value = std::is_same<function_prototype<F>, function_prototype<R>>::value;
};

template <typename F, typename R>
constexpr bool same_prototype_v = same_prototype<F, R>::value;
