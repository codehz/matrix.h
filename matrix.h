#pragma once
#include "matrix_traits.h"
#include "util.h"
#include <iomanip>
#include <iostream>
#include <memory>
#include <utility>

namespace sm {

template <typename T, typename Item, typename Line>
struct matrix_printer {
    const T &t;
    Item item;
    Line line;

    friend std::ostream &operator<<(std::ostream &os, const matrix_printer &printer)
    {
        printer.t.for_each_row([&os, item = printer.item, line = printer.line ](
            row_helper<typename T::number_t, T::rows, T::cols, typename T::store_t, true> r,
            std::size_t row) {
            r.for_each([&os, item, row](typename T::number_t num, std::size_t col) {
                item(os, num, row, col);
            });
            line(os, row);
        });
        return os;
    }
};

template <typename Number, std::size_t Rows, std::size_t Cols,
          template <typename, std::size_t, std::size_t,
                    template <typename, std::size_t, std::size_t> typename> typename RealType,
          template <typename, std::size_t, std::size_t> typename StoreType>
struct matrix_base : public matrix_traits<Number, Rows, Cols> {
    template <typename Number_, std::size_t Rows_, std::size_t Cols_,
              template <typename, std::size_t, std::size_t> typename StoreType_>
    using real_template = RealType<Number_, Rows_, Cols_, StoreType_>;
    using real_t = RealType<Number, Rows, Cols, StoreType>;
    using store_t = StoreType<Number, Rows, Cols>;
    template <typename Number_, std::size_t Rows_, std::size_t Cols_>
    using store_type = StoreType<Number_, Rows_, Cols_>;
    using real_array_t = RealType<Number, Rows, Cols, store_array>;

    static constexpr range<0, Rows> allRows = {};
    static constexpr range<0, Cols> allCols = {};

    store_t data;

    matrix_base() noexcept {}

    template <typename T,
              typename = typename std::enable_if<T::rows == Rows && T::cols == Cols>::type>
    matrix_base(const T &t)
    {
        for_each_with([](Number source) { return source; }, t);
    }

    matrix_base(const Number (&arr)[Rows][Cols])
        : matrix_base(array_wrapper<const Number, Rows, Cols>(arr))
    {
    }

    template <typename T, typename = typename std::enable_if<
                              std::is_constructible<store_t, const T &>::value>::type>
    matrix_base(matrix_no_copy_t no_copy, const T &proxy) : data(proxy)
    {
    }

    row_helper<Number, Rows, Cols, store_t, true> operator[](std::size_t row) const noexcept
    {
        return {data, row, direct_access<Number, store_t, true>};
    }

    row_helper<Number, Rows, Cols, store_t, false> operator[](std::size_t row) noexcept
    {
        return {data, row, direct_access<Number, store_t, false>};
    }

    template <typename T>
    std::enable_if_t<T::rows == Rows && T::cols == Cols, real_t &> set(const T &t)
    {
        for_each_with([](Number &target, Number source) { target = source; }, t);
        return static_cast<real_t &>(*this);
    }

    template <typename T, typename EX = std::decay_t<T>>
    std::enable_if_t<EX::rows == Rows && EX::cols == Cols> swap(T &&t)
    {
        auto temp = t.clone_to_array();
        t.set(*this);
        set(temp);
    }

    template <typename F>
    real_t &for_each_row_1(F f)
    {
        static_for<std::size_t, 0, Rows>([&, this](auto row) { f((*this)[row]); });
        return static_cast<real_t &>(*this);
    }

    template <typename F>
    const real_t &for_each_row_1(F f) const
    {
        static_for<std::size_t, 0, Rows>([&, this](auto row) { f((*this)[row]); });
        return static_cast<const real_t &>(*this);
    }

    template <typename F>
    real_t &for_each_row_2(F f)
    {
        static_for<std::size_t, 0, Rows>([&, this](auto row) { f((*this)[row], row); });
        return static_cast<real_t &>(*this);
    }

    template <typename F>
    const real_t &for_each_row_2(F f) const
    {
        static_for<std::size_t, 0, Rows>([&, this](auto row) { f((*this)[row], row); });
        return static_cast<const real_t &>(*this);
    }

    template <typename F>
    real_t &for_each_row(F f)
    {
        static_for<std::size_t, 0, Rows>([&, this](auto row) {
            static_if(function_traits<F>::arity == 2) f((*this)[row], row);
            else static_if(function_traits<F>::arity == 1) f((*this)[row]);
            else throw "for_each_with failed";
        });
        return static_cast<real_t &>(*this);
    }

    template <typename F>
    const real_t &for_each_row(F f) const
    {
        static_for<std::size_t, 0, Rows>([&, this](auto row) {
            static_if(function_traits<F>::arity == 2) f((*this)[row], row);
            else static_if(function_traits<F>::arity == 1) f((*this)[row]);
            else throw "for_each_row(const) failed";
        });
        return static_cast<const real_t &>(*this);
    }

    template <typename F>
    real_t &for_each(F f)
    {
        static_for<std::size_t, 0, Rows>([&, this](auto row) {
            static_for<std::size_t, 0, Cols>([&, this](auto col) {
                static_if(std::is_same<typename function_traits<F>::return_type, void>::value)
                {
                    static_if(function_traits<F>::arity == 1) f(data[row][col]);
                    else static_if(function_traits<F>::arity == 0) f();
                    else static_if(function_traits<F>::arity == 3) f(data[row][col], row, col);
                    else throw "for_each failed";
                }
                else
                {
                    static_if(function_traits<F>::arity == 0) data[row][col] = f();
                    else static_if(function_traits<F>::arity == 1) data[row][col] =
                        f(data[row][col]);
                    else static_if(function_traits<F>::arity == 2) data[row][col] = f(row, col);
                    else static_if(function_traits<F>::arity == 3) data[row][col] =
                        f(data[row][col], row, col);
                    else throw "for_each failed";
                }
            });
        });
        return static_cast<real_t &>(*this);
    }

    template <typename F>
    real_t &for_each(F f) const
    {
        static_for<std::size_t, 0, Rows>([&, this](auto row) {
            static_for<std::size_t, 0, Cols>([&, this](auto col) {
                static_if(function_traits<F>::arity == 1) f(data[row][col]);
                else static_if(function_traits<F>::arity == 0) f();
                else static_if(function_traits<F>::arity == 3) f(data[row][col], row, col);
                else throw "for_each(const) failed";
            });
        });
        return static_cast<real_t &>(*this);
    }

    template <typename F, typename... R>
    real_t &for_each_with(F f, R... r)
    {
        static_assert(((R::rows == Rows && R::cols == Cols) && ...), "");
        static_for<std::size_t, 0, Rows>([&, this](auto row) {
            static_for<std::size_t, 0, Cols>([&, this](auto col) {
                static_if(
                    same_prototype_v<F, Number(std::size_t, std::size_t, typename R::number_t...)>)
                    data[row][col] = f(row, col, r[row][col]...);
                else static_if(same_prototype_v<F, Number(typename R::number_t...)>)
                    data[row][col] = f(r[row][col]...);
                else static_if(same_prototype_v<F, Number(Number, std::size_t, std::size_t,
                                                          typename R::number_t...)>)
                    data[row][col] = f(data[row][col], row, col, r[row][col]...);
                else static_if(same_prototype_v<F, Number(Number, typename R::number_t...)>)
                    data[row][col] = f(data[row][col], r[row][col]...);
                else static_if(same_prototype_v<F, void(Number &, std::size_t, std::size_t,
                                                        typename R::number_t...)> ||
                               same_prototype_v<F, void(Number, std::size_t, std::size_t,
                                                        typename R::number_t...)>)
                    f(data[row][col], row, col, r[row][col]...);
                else static_if(same_prototype_v<F, void(Number &, typename R::number_t...)> ||
                               same_prototype_v<F, void(Number, typename R::number_t...)>)
                    f(data[row][col], r[row][col]...);
                else throw "for_each_with failed";
            });
        });
        return static_cast<real_t &>(*this);
    }

    template <typename F, typename... R>
    real_t &for_each_with(F f, R... r) const
    {
        static_assert(((R::rows == Rows && R::cols == Cols) && ...), "");
        static_for<std::size_t, 0, Rows>([&, this](auto row) {
            static_for<std::size_t, 0, Cols>([&, this](auto col) {
                static_if(same_prototype_v<F, void(Number, std::size_t, std::size_t,
                                                   typename R::number_t...)>)
                    f(data[row][col], row, col, r[row][col]...);
                else static_if(same_prototype_v<F, void(Number, typename R::number_t...)>)
                    f(data[row][col], r[row][col]...);
                else throw "for_each_with(const) failed";
            });
        });
        return static_cast<real_t &>(*this);
    }

    auto mirror()
    {
        return make_proxy<Number, Cols, Rows, RealType, matrix_base, false>(
            *this, mirror_access<Number, matrix_base, true>,
            transpose_access<Number, matrix_base, false>);
    }

    auto mirror() const
    {
        return make_proxy<Number, Cols, Rows, RealType, matrix_base, true>(
            *this, mirror_access<Number, matrix_base, true>,
            transpose_access<Number, matrix_base, true>);
    }

    auto transpose()
    {
        return make_proxy<Number, Cols, Rows, RealType, matrix_base, false>(
            *this, transpose_access<Number, matrix_base, true>,
            transpose_access<Number, matrix_base, false>);
    }

    auto transpose() const
    {
        return make_proxy<Number, Cols, Rows, RealType, matrix_base, true>(
            *this, transpose_access<Number, matrix_base, true>,
            transpose_access<Number, matrix_base, true>);
    }

    template<std::size_t nRows, std::size_t nCols>
    auto vexpand(const matrix_size<nRows, nCols> &) {
        return make_proxy<Number, nRows, nCols, RealType, matrix_base, false>(
            *this, expand_access<Number, matrix_base, true>,
            expand_access<Number, matrix_base, false>);
    }

    template<std::size_t nRows, std::size_t nCols>
    auto vexpand(const matrix_size<nRows, nCols> &) const {
        return make_proxy<Number, nRows, nCols, RealType, matrix_base, true>(
            *this, expand_access<Number, matrix_base, true>,
            expand_access<Number, matrix_base, true>);
    }

    template <typename T>
    auto operator+(const T &rhs) const
    {
        static_assert(T::rows == Rows && T::cols == Cols, "");
        return real_array_t{}.for_each_with([](Number a, Number b) { return a + b; }, *this, rhs);
    }

    template <typename T>
    auto operator-(const T &rhs) const
    {
        static_assert(T::rows == Rows && T::cols == Cols, "");
        return real_array_t{}.for_each_with([](Number a, Number b) { return a - b; }, *this, rhs);
    }

    template <typename T>
    real_t &operator+=(const T &rhs)
    {
        static_assert(T::rows == Rows && T::cols == Cols, "");
        return for_each_with([](Number target, Number b) { return target + b; }, rhs);
    }

    template <typename T>
    real_t &operator-=(const T &rhs)
    {
        static_assert(T::rows == Rows && T::cols == Cols, "");
        return for_each_with([](Number target, Number b) { return target - b; }, rhs);
    }

    template <std::size_t N, template <typename, std::size_t, std::size_t> typename RStoreType>
    auto operator*(const RealType<Number, Cols, N, RStoreType> &rhs) const
    {
        return RealType<Number, Rows, N, store_array>{}.for_each(
            [this, &rhs](Number &target, std::size_t row, std::size_t col) {
                target = 0;
                static_for<std::size_t, 0, Cols>(
                    [&, this](auto n) { target += (*this)[row][n] * rhs[n][col]; });
            });
    }

    template <template <typename, std::size_t, std::size_t> typename RStoreType>
    real_t &operator*=(const RealType<Number, Cols, Cols, RStoreType> &rhs)
    {
        return for_each_with([](Number source) { return source; }, (*this) * rhs);
    }

    real_array_t operator*(Number num) const
    {
        return transform([num](Number src) { return src * num; });
    }

    real_array_t operator/(Number num) const
    {
        return transform([num](Number src) { return src / num; });
    }

    real_t &operator*=(Number num)
    {
        return for_each([num](Number src) { return src * num; });
    }

    real_t &operator/=(Number num)
    {
        return for_each([num](Number src) { return src / num; });
    }

    friend real_array_t operator*(Number num, const matrix_base &m) { return m * num; }

    template <typename F>
    std::enable_if_t<same_prototype_v<F, Number(Number, std::size_t, std::size_t)>, real_array_t>
    transform(F f) const
    {
        return real_array_t{}.for_each_with([f](Number &num, std::size_t row, std::size_t col,
                                                Number src) { num = f(src, row, col); },
                                            *this);
    }

    template <typename F>
    std::enable_if_t<same_prototype_v<F, Number(Number)>, real_array_t> transform(F f) const
    {
        return transform([f](Number num, std::size_t, std::size_t) { return f(num); });
    }

    Number norm()
    {
        Number ret;
        for_each([&ret](Number num) { ret += num * num; });
        return ret;
    }

    auto &normalize() { return (*this) /= norm(); }

    template <typename RowRange, typename ColRange>
    auto operator()(const RowRange &rr, const ColRange &cc)
    {
        static_assert(RowRange::max <= Rows && ColRange::max <= Cols, "");
        return make_proxy<Number, RowRange::length, ColRange::length, RealType,
                          StoreType<Number, Rows, Cols>, false>(
            this->data,
            range_transform<RowRange, ColRange>::template func<Number,
                                                               StoreType<Number, Rows, Cols>, true>,
            range_transform<RowRange,
                            ColRange>::template func<Number, StoreType<Number, Rows, Cols>, false>);
    }

    template <typename RowRange, typename ColRange>
    auto operator()(const RowRange &rr, const ColRange &cc) const
    {
        static_assert(RowRange::max <= Rows && ColRange::max <= Cols, "");
        return make_proxy<Number, RowRange::length, ColRange::length, RealType,
                          StoreType<Number, Rows, Cols>, true>(
            this->data,
            range_transform<RowRange, ColRange>::template func<Number,
                                                               StoreType<Number, Rows, Cols>, true>,
            range_transform<RowRange,
                            ColRange>::template func<Number, StoreType<Number, Rows, Cols>, true>);
    }

    template <template <typename, std::size_t, std::size_t> typename TargetStoreType = StoreType>
    RealType<Number, Rows, Cols, TargetStoreType> clone() const
    {
        return RealType<Number, Rows, Cols, TargetStoreType>{static_cast<const real_t &>(*this)};
    }

    RealType<Number, Rows, Cols, store_array> clone_to_array() const
    {
        return RealType<Number, Rows, Cols, store_array>{static_cast<const real_t &>(*this)};
    }

    template <typename Item, typename Line>
    std::enable_if_t<
        same_prototype_v<Item, void(std::ostream &, Number, std::size_t, std::size_t)> &&
            same_prototype_v<Line, void(std::ostream &, std::size_t)>,
        matrix_printer<real_t, Item, Line>>
    printer(Item item, Line line) const noexcept
    {
        return {static_cast<const real_t &>(*this), item, line};
    };

    friend std::ostream &operator<<(std::ostream &os, const matrix_base &m)
    {
        m.for_each_row([&os](row_helper<Number, Rows, Cols, store_t, true> r) {
            r.for_each([&os](Number num) { os << num << '\t'; });
            os << std::endl;
        });
        return os;
    }

    void _id_matrix_base() noexcept {}
};

template <typename Number, std::size_t Rows, std::size_t Cols,
          template <typename, std::size_t, std::size_t> typename StoreType = sm::store_array>
struct matrix : public matrix_base<Number, Rows, Cols, matrix, StoreType> {
    using super_t = typename find_type<decltype(&matrix::_id_matrix_base)>::type;

    template <typename... T>
    matrix(T... t) : super_t(std::forward<T>(t)...)
    {
    }
    matrix(const Number (&arr)[super_t::rows][super_t::cols])
        : super_t(std::forward<const Number (&)[super_t::rows][super_t::cols]>(arr))
    {
    }
};

template <typename Matrix, bool hori = Matrix::rows == 1>
struct vector_base : Matrix {
    static constexpr std::size_t length = hori ? Matrix::cols : Matrix::rows;
    using super_t = Matrix;

    template <typename... T>
    vector_base(T... t) : super_t(std::forward<T>(t)...)
    {
    }

    vector_base(const typename Matrix::number_t (&arr)[super_t::rows][super_t::cols])
        : super_t(
              std::forward<const typename Matrix::number_t (&)[super_t::rows][super_t::cols]>(arr))
    {
    }

    auto &operator()(std::size_t pos)
    {
        static_if(hori) return (*this)[0][pos];
        else return (*this)[pos][0];
    }

    auto operator()(std::size_t pos) const
    {
        static_if(hori) return (*this)[0][pos];
        else return (*this)[pos][0];
    }

    template <typename Range>
    auto operator()(const Range &r)
    {
        static_if(hori) return this->Matrix::operator()(this->allRows, r);
        else return this->Matrix::operator()(r, this->allCols);
    }

    template <typename Range>
    auto operator()(const Range &r) const
    {
        static_if(hori) return this->Matrix::operator()(this->allRows, r);
        else return this->Matrix::operator()(r, this->allCols);
    }

    template <typename RowRange, typename ColRange>
    auto operator()(const RowRange &rr, const ColRange &cc)
    {
        return this->Matrix::operator()(rr, cc);
    }

    template <typename RowRange, typename ColRange>
    auto operator()(const RowRange &rr, const ColRange &cc) const
    {
        return this->Matrix::operator()(rr, cc);
    }

    void _id_vector_base() noexcept {}
};

template <typename Number, std::size_t Length,
          template <typename, std::size_t, std::size_t> typename StoreType>
struct matrix<Number, 1, Length, StoreType>
    : public vector_base<matrix_base<Number, 1, Length, matrix, StoreType>> {
    using super_t = typename find_type<decltype(&matrix::_id_vector_base)>::type;

    template <typename... T>
    matrix(T... t) : super_t(std::forward<T>(t)...)
    {
    }
    matrix(const Number (&arr)[super_t::rows][super_t::cols])
        : super_t(std::forward<const Number (&)[super_t::rows][super_t::cols]>(arr))
    {
    }
};

template <typename Number, std::size_t Length,
          template <typename, std::size_t, std::size_t> typename StoreType>
struct matrix<Number, Length, 1, StoreType>
    : public vector_base<matrix_base<Number, Length, 1, matrix, StoreType>> {
    using super_t = typename find_type<decltype(&matrix::_id_vector_base)>::type;

    template <typename... T>
    matrix(T... t) : super_t(std::forward<T>(t)...)
    {
    }
    matrix(const Number (&arr)[super_t::rows][super_t::cols])
        : super_t(std::forward<const Number (&)[super_t::rows][super_t::cols]>(arr))
    {
    }
};

template <typename Number, template <typename, std::size_t, std::size_t> typename StoreType>
struct matrix<Number, 1, 1, StoreType>
    : public vector_base<matrix_base<Number, 1, 1, matrix, StoreType>> {
    using super_t = typename find_type<decltype(&matrix::_id_vector_base)>::type;

    template <typename... T>
    matrix(T... t) : super_t(std::forward<T>(t)...)
    {
    }
    matrix(const Number (&arr)[super_t::rows][super_t::cols])
        : super_t(std::forward<const Number (&)[super_t::rows][super_t::cols]>(arr))
    {
    }

    matrix(const Number &num)
        : super_t({{num}})
    {
    }

    auto &operator()(std::size_t pos) { return super_t::operator()(pos); }

    auto operator()(std::size_t pos) const { return super_t::operator()(pos); }

    template <typename Range>
    auto operator()(const Range &r)
    {
        return super_t::operator()(r);
    }

    template <typename Range>
    auto operator()(const Range &r) const
    {
        return super_t::operator()(r);
    }

    template <typename RowRange, typename ColRange>
    auto operator()(const RowRange &rr, const ColRange &cc)
    {
        return this->Matrix::operator()(rr, cc);
    }

    template <typename RowRange, typename ColRange>
    auto operator()(const RowRange &rr, const ColRange &cc) const
    {
        return this->Matrix::operator()(rr, cc);
    }

    auto &operator()() { return (*this)[0][0]; }

    auto operator()() const { return (*this)[0][0]; }
};

template <typename Number, std::size_t Length,
          template <typename, std::size_t, std::size_t> typename StoreType>
struct matrix<Number, Length, Length, StoreType>
    : public matrix_base<Number, Length, Length, matrix, StoreType> {
    using super_t = typename find_type<decltype(&matrix::_id_matrix_base)>::type;

    template <typename... T>
    matrix(T... t) : super_t(std::forward<T>(t)...)
    {
    }
    matrix(const Number (&arr)[super_t::rows][super_t::cols])
        : super_t(std::forward<const Number (&)[super_t::rows][super_t::cols]>(arr))
    {
    }

    auto LUPDecomposition()
    {
        typename super_t::real_array_t l, u, p;
        auto &a = *this;
        auto zero = [](Number &num) { num = 0; };
        auto eye = [](std::size_t row, std::size_t col) -> Number { return row == col; };
        auto gabs = [](Number num) { return num >= 0 ? num : -num; };
        l.for_each(eye);
        u.for_each(zero);
        p.for_each(eye);
        static_for<std::size_t, 0, Length>([&](auto i) {
            auto max_j = i.value;
            static_for<std::size_t, i, Length>([&, a = a.transform(gabs) ](auto j) {
                if (a[j][i] > a[max_j][i]) max_j = j;
            });
            if (max_j != i)
                static_for<std::size_t, 0, Length>(
                    [&](auto k) { std::swap(p[i][k], p[max_j][k]); });
        });

        auto aprime = (p * a);

        static_for<std::size_t, 0, Length>([&](auto i) {
            static_for<std::size_t, 0, Length>([&](auto j) {
                if (j <= i) {
                    Number s = 0;
                    static_for<std::size_t, 0, j>([&](auto k) { s += l[j][k] * u[k][i]; });
                    u[j][i] = aprime[j][i] - s;
                }
                else if (j >= i) {
                    Number s = 0;
                    static_for<std::size_t, 0, i>([&](auto k) { s += l[j][k] * u[k][i]; });
                    l[j][i] = (aprime[j][i] - s) / u[i][i];
                }
            });
        });

        return std::make_tuple(l, u, p);
    }

    template <template <typename, std::size_t, std::size_t> typename RStoreType>
    auto solve_equaltions(const typename super_t::template real_template<
                          typename super_t::number_t, super_t::rows, 1, RStoreType> &b)
    {
        typename super_t::real_array_t l, u, p;
        std::tie(l, u, p) = LUPDecomposition();
        typename super_t::template real_template<typename super_t::number_t, super_t::rows, 1,
                                                 StoreType>
            x, y;
        static_for<std::size_t, 0, Length>([&, tb = p * b ](auto I) {
            constexpr auto i = I.value;
            static_if(i == 0) y(i) = tb(i);
            else y(i) = tb(i) - (l(sm::range_one_v<i>, sm::range_v<0, i>) * y(sm::range_v<0, i>))();
        });
        static_for<std::size_t, 0, Length, true>([&](auto I) {
            constexpr auto i = I.value;
            static_if(i + 1 == Length) x(i) = y(i) / u[i][i];
            else x(i) = (y(i) - (u(sm::range_one_v<i>, sm::range_v<i + 1, Length>) *
                                 x(sm::range_v<i + 1, Length>))()) /
                        u[i][i];
        });
        return x;
    }
};

template <typename Number, std::size_t Length,
          template <typename, std::size_t, std::size_t> typename StoreType = store_array>
using vector = matrix<Number, 1, Length, StoreType>;

template <typename Number, std::size_t Length,
          template <typename, std::size_t, std::size_t> typename StoreType = store_array>
using vector_vert = matrix<Number, Length, 1, StoreType>;

template <typename Number, 
          template <typename, std::size_t, std::size_t> typename StoreType = store_array>
using singularity = matrix<Number, 1, 1, StoreType>;
}
