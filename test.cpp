//#define NTYPE
#include "dynamic_template.h"
#include "matrix.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>
#ifndef NTYPE
#include <cxxabi.h>
template <class T>
void print_type_of(T const &o)
{
    char const *mangled = typeid(o).name();
    int status;
    char *demangled = abi::__cxa_demangle(mangled, 0, 0, &status);
    std::cout << "Mangled: " << mangled << "\n";
    if (status) {
        std::cout << "Demangling failed with status " << status << "\n";
    }
    else {
        std::cout << "Demangled: " << demangled << "\n";
    }
    std::cout << std::endl;
}
#else
template <typename T>
void print_type_of(T const &o)
{
}
#endif

int main()
{
    std::size_t n = 0;
    std::cout << "n: ";
    std::cin >> n;
    dynamic_template<std::size_t, 2, 4>::eval(
        n,
        [](auto N) {
            sm::matrix<double, N, N> a{};
            sm::vector_vert<double, N> b{};
            print_type_of(a);
            print_type_of(b);

            std::cout << "a:" << std::endl;
            static_for<std::size_t, 0, N>([&](auto i) {
                static_for<std::size_t, 0, N>([&](auto j) { std::cin >> a[i][j]; });
            });
            std::cout << "b:" << std::endl;
            static_for<std::size_t, 0, N>([&](auto i) { std::cin >> b(i.value); });

            auto item = [](std::ostream &os, double num, std::size_t, std::size_t) {
                os << std::setw(12) << std::left << num;
            };

            auto line = [](std::ostream &os, std::size_t) { os << std::endl; };

            std::cout << std::setprecision(6) << std::showpos;

            std::cout << "a: \n" << a.printer(item, line) << std::endl;
            std::cout << "b: \n" << b.printer(item, line) << std::endl;
            std::cout << "solve: \n" << a.solve_equaltions(b);
        });
}
