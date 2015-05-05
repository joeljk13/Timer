#define BOOST
#include "timer.hpp"

#include <chrono>

#include <algorithm>
#include <iostream>
#include <deque>
#include <vector>
#include <random>

#include <ctime>
#include <cstddef>
#include <cstdlib>

namespace {

struct func1 {
    template <class T>
    __attribute__ ((__always_inline__,
                    __hot__))
    void operator ()(T t, std::size_t bit, T value)
    {
        volatile T tmp = value == static_cast<T>(1u) ? 
            t & ~(static_cast<T>(1u) << bit) :
            t | static_cast<T>(1u) << bit;
    }
};

struct func2 {
    template <class T>
    __attribute__ ((__always_inline__,
                    __hot__))
    void operator ()(T t, std::size_t bit, T value)
    {
        volatile T tmp = t & ~(static_cast<T>(1u) << bit) | value << bit;
    }
};

}

int main()
{
    std::chrono::seconds duration (5);
    timer::default_alpha = 0.1l;
    std::srand(std::time(NULL));
    unsigned int r, v;
    std::size_t b;
#define SET do { \
    r = (unsigned int)std::rand(); \
    b = (std::size_t)std::rand() % 32; \
    v = (unsigned int)std::rand() % 2; \
    std::cout \
        << "r = " << r \
        << "\nb = " << b \
        << "\nv = " << v \
        << "\n!!(r & 1u << b) == v = " << (!!(r & 1u << b) == v) \
        << std::endl; \
} while (0)
    SET;
    timer branch (func1(), duration, r, b, v);
    std::cout
        << "Branch ratio: "
        << branch.get_ratio()
        << std::endl;
    SET;
    branch.measure(func1(), duration, r, b, v);
    std::cout
        << "Branch ratio: "
        << branch.get_ratio()
        << std::endl;
    SET;
    branch.measure(func1(), duration, r, b, v);
    std::cout
        << "Branch ratio: "
        << branch.get_ratio()
        << std::endl;
    SET;
    branch.measure(func1(), duration, r, b, v);
    std::cout
        << "Branch ratio: "
        << branch.get_ratio()
        << std::endl;
    SET;
    timer bits (func2(), duration, r, b, v);
    std::cout
        << "Bits ratio: "
        << bits.get_ratio()
        << std::endl;
    SET;
    bits.measure(func2(), duration, r, b, v);
    std::cout
        << "Bits ratio: "
        << bits.get_ratio()
        << std::endl;
    SET;
    bits.measure(func2(), duration, r, b, v);
    std::cout
        << "Bits ratio: "
        << bits.get_ratio()
        << std::endl;
    SET;
    bits.measure(func2(), duration, r, b, v);
    std::cout
        << "Bits ratio: "
        << bits.get_ratio()
        << std::endl;
    std::cout
        << "\nAssign by branch: "
        << branch.get_ratio()
        << "\nAssign by bits: "
        << bits.get_ratio()
        // << timer::compare(func1(), func2(), duration, r, b, v)
        << std::endl;
    return 0;
}

/*
int main()
{
    const int bins = 50;
    const int count = 1000;
    const int repetitions = 1000;
    const char marker = '*';
    int **array = new int *[bins];
    for (int i = 0; i < bins; ++i) {
        array[i] = new int[count];
    }
    dist(&fn, count, array, bins, repetitions);
    std::cout << std::endl << "-----Actual Output-----" << std::endl <<
        std::endl;
    for (int i = 0; i < bins; ++i) {
        for (int j = 0; j < count; ++j) {
            if (array[i][j] > 0) {
                std::cout << marker;
            }
        }
        std::cout << std::endl << std::endl;
    }
    for (int i = 0; i < bins; ++i) {
        delete [] array[i];
    }
    delete [] array;
    return 0;
}
*/
