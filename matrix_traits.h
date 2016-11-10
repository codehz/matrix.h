#pragma once
#include "dynamic_template.h"
#include "util.h"
#include <cassert>

namespace sm {

struct matrix_no_copy_t {
};

constexpr matrix_no_copy_t matrix_no_copy{};

template <std::size_t Rows, std::size_t Cols>
struct matrix_size {
    static constexpr std::size_t rows = Rows;
    static constexpr std::size_t cols = Cols;
};

template <std::size_t Rows, std::size_t Cols>
constexpr matrix_size<Rows, Cols> matrix_size_v{};

template <typename Number, std::size_t Rows, std::size_t Cols>
struct matrix_traits : matrix_size<Rows, Cols> {
    using number_t = Number;
    using number_ref = Number &;
};

template <typename T>
struct matrix_info : T {
};

template <typename Number, std::size_t Rows, std::size_t Cols>
struct matrix_info<Number[Rows][Cols]> : matrix_traits<Number, Rows, Cols> {
};

template <typename Number, typename SrcStoreType, bool is_const>
using access_fn = const_switch_t<is_const, Number> (*)(const_switch_t<is_const, SrcStoreType &> Src,
                                                       std::size_t row, std::size_t col);

template <typename Number, typename SrcStoreType, bool is_const = false>
const_switch_t<is_const, Number> direct_access(const_switch_t<is_const, SrcStoreType &> Src,
                                               std::size_t row, std::size_t col)
{
    return Src[row][col];
}

template <typename Number, typename SrcStoreType, bool is_const = false>
const_switch_t<is_const, Number> transpose_access(const_switch_t<is_const, SrcStoreType &> Src,
                                                  std::size_t row, std::size_t col)
{
    return Src[col][row];
}

template <typename Number, typename SrcStoreType, bool is_const = false>
const_switch_t<is_const, Number> expand_access(const_switch_t<is_const, SrcStoreType &> Src,
                                               std::size_t row, std::size_t col)
{
    return Src[row % matrix_info<SrcStoreType>::rows][col % matrix_info<SrcStoreType>::cols];
}

template <typename Number, typename SrcStoreType, bool is_const = false>
const_switch_t<is_const, Number> mirror_access(const_switch_t<is_const, SrcStoreType &> Src,
                                               std::size_t row, std::size_t col)
{
    return Src[matrix_info<SrcStoreType>::rows - row - 1]
              [matrix_info<SrcStoreType>::cols - col - 1];
}

template <typename Number, std::size_t Rows, std::size_t Cols, typename SrcStoreType,
          bool is_const = false>
struct row_helper : public matrix_traits<Number, Rows, Cols> {
    using super_t = matrix_traits<Number, Rows, Cols>;
    using store_t = SrcStoreType;
    using data_t = const_switch_t<is_const, store_t &>;
    data_t data;
    std::size_t current_row;
    access_fn<Number, store_t, is_const> access;

    row_helper(data_t data, std::size_t row,
               access_fn<Number, store_t, is_const> f = &direct_access<Number, store_t, is_const>)
        : data(data), current_row(row), access(f)
    {
    }

    Number operator[](std::size_t col) const { return access(data, current_row, col); }

    template <bool cond = !is_const>
    std::enable_if_t<cond, Number &> operator[](std::size_t col)
    {
        return access(data, current_row, col);
    }

    template <typename F>
    std::enable_if_t<same_prototype_v<F, void(Number, std::size_t)>> for_each(F f) const
    {
        static_for<std::size_t, 0, Cols>([&](auto col) { f(data[current_row][col], col); });
    }

    template <typename F, bool cond = !is_const>
    std::enable_if_t<cond && same_prototype_v<F, void(Number &, std::size_t)>> for_each(F f) const
    {
        static_for<std::size_t, 0, Cols>([&](auto col) { f(data[current_row][col], col); });
    }

    template <typename F>
    std::enable_if_t<same_prototype_v<F, void(Number)>> for_each(F f) const
    {
        for_each([f](Number num, std::size_t col) { return f(num); });
    }

    template <typename F>
    std::enable_if_t<same_prototype_v<F, void(Number &)>> for_each(F f)
    {
        for_each([f](Number &num, std::size_t col) { return f(num); });
    }
};

template <typename Number, std::size_t Rows, std::size_t Cols>
using store_array = Number[Rows][Cols];

template <typename Number, std::size_t Rows, std::size_t Cols>
struct array_wrapper : public matrix_traits<Number, Rows, Cols> {
    const Number (&source)[Rows][Cols];

    array_wrapper(const Number (&source)[Rows][Cols]) : source(source) {}

    row_helper<Number, Rows, Cols, store_array<Number, Rows, Cols>, true>
    operator[](std::size_t row) const
    {
        return {source, row};
    }
};

template <typename Number, std::size_t Rows, std::size_t Cols, typename SrcStoreType, bool is_const>
struct store_proxy : matrix_traits<Number, Rows, Cols> {
    using store_t = const_switch_t<is_const, SrcStoreType &>;
    store_t ref;
    const access_fn<Number, store_t, true> const_access;
    const access_fn<Number, store_t, is_const> access;

    template <typename T>
    std::enable_if_t<T::rows == Rows && T::cols == Cols, store_proxy &> set(T t)
    {
        static_for<std::size_t, 0, Rows>([&, this](auto row) {
            static_for<std::size_t, 0, Cols>([&, this](auto col) { ref[row][col] = t[row][col]; });
        });
    }

    template <bool cond = is_const, typename = std::enable_if_t<cond>>
    store_proxy(store_t ref, access_fn<Number, store_t, true> f)
        : ref(ref), const_access(f), access(f)
    {
    }

    template <bool cond = !is_const, typename = std::enable_if_t<cond>>
    store_proxy(store_t ref, access_fn<Number, store_t, true> const_f,
                access_fn<Number, store_t, false> f)
        : ref(ref), const_access(const_f), access(f)
    {
    }

    row_helper<Number, Rows, Cols, SrcStoreType, true> operator[](std::size_t row) const
    {
        return {ref, row, const_access};
    }

    template <bool cond = !is_const>
    std::enable_if_t<cond, row_helper<Number, Rows, Cols, SrcStoreType, false>>
    operator[](std::size_t row)
    {
        return {ref, row, access};
    }
};

template <typename SrcStoreType, bool is_const>
struct store_proxy_t {
    template <typename Number, std::size_t Rows, std::size_t Cols>
    using type = store_proxy<Number, Rows, Cols, SrcStoreType, is_const>;
};

template <typename Number, std::size_t Rows, std::size_t Cols,
          template <typename, std::size_t, std::size_t,
                    template <typename, std::size_t, std::size_t> typename> typename RealType,
          typename SrcStoreType, bool is_const>
RealType<Number, Rows, Cols, store_proxy_t<SrcStoreType, is_const>::template type>
make_proxy(const_switch_t<is_const, SrcStoreType &> src,
           access_fn<Number, SrcStoreType, true> const_f,
           access_fn<Number, SrcStoreType, is_const> f)
{
    return {matrix_no_copy,
            store_proxy<Number, Rows, Cols, SrcStoreType, is_const>(src, const_f, f)};
};

template <typename... RangesType>
struct range_pack;

template <typename... RangesType>
constexpr range_pack<RangesType...> range_pack_v{};

template <std::size_t Low, std::size_t High>
struct range {
    static_assert(Low < High, "Low should < High");
    static constexpr std::size_t low = Low;
    static constexpr std::size_t hight = High;

    static constexpr std::size_t length = High - Low;
    static constexpr std::size_t start = 0;
    static constexpr std::size_t max = High;

    static std::size_t transform(std::size_t pos)
    {
        assert(pos < length);
        return Low + pos;
    }

    template <std::size_t nLows, std::size_t nHighs>
    constexpr auto operator+(range<nLows, nHighs>) const
    {
        return range_pack_v<range<nLows, nHighs>, range<Low, High>>;
    }

    template <std::size_t... nrLows, std::size_t... nrHighs,
              template <std::size_t, std::size_t> typename... nrT>
    constexpr auto operator+(nrT<nrLows, nrHighs>...) const
    {
        return range_pack_v<range<nrLows, nrHighs>..., range<Low, High>>;
    }
};

template <std::size_t Low, std::size_t High>
constexpr range<Low, High> range_v{};

template <std::size_t One>
struct range_one : range<One, One + 1> {
};

template <std::size_t One>
constexpr range_one<One> range_one_v{};

template <std::size_t Lows, std::size_t Highs, template <std::size_t, std::size_t> typename T>
struct range_pack<T<Lows, Highs>> : range<Lows, Highs> {
};

template <std::size_t fLows, std::size_t... rLows, std::size_t fHighs, std::size_t... rHighs,
          template <std::size_t, std::size_t> typename fT,
          template <std::size_t, std::size_t> typename... rT>
struct range_pack<fT<fLows, fHighs>, rT<rLows, rHighs>...> : range_pack<rT<rLows, rHighs>...> {
    using super_t = range_pack<rT<rLows, rHighs>...>;
    static_assert(fLows > super_t::max, "please sort");
    static constexpr std::size_t length = super_t::length + fHighs - fLows;
    static constexpr std::size_t start = super_t::start + super_t::length;
    static constexpr std::size_t max = fHighs;

    static std::size_t transform(std::size_t pos)
    {
        assert(pos < length);
        if (pos < start) return super_t::transform(pos);
        return pos - start + fLows;
    }

    template <std::size_t nLows, std::size_t nHighs>
    constexpr auto operator+(range<nLows, nHighs>) const
    {
        return range_pack_v<range<nLows, nHighs>, range<fLows, fHighs>, rT<rLows, rHighs>...>;
    }

    template <std::size_t... nrLows, std::size_t... nrHighs,
              template <std::size_t, std::size_t> typename... nrT>
    constexpr auto operator+(nrT<nrLows, nrHighs>...) const
    {
        return range_pack_v<range<nrLows, nrHighs>..., range<fLows, fHighs>, rT<rLows, rHighs>...>;
    }
};

template <typename RowRange, typename ColRange>
struct range_transform {
    template <typename Number, typename SrcStoreType, bool is_const>
    static const_switch_t<is_const, Number> func(const_switch_t<is_const, SrcStoreType &> Src,
                                                 std::size_t row, std::size_t col)
    {
        return Src[RowRange::transform(row)][ColRange::transform(col)];
    }
};
}
