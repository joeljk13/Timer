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

typedef unsigned int T;

T func1(T t, std::size_t bit) {
    return t >> bit & 1u;
}

T func2(T t, std::size_t bit) {
    return t & 1u << bit ? 1 : 0;
}

T func3(T t, std::size_t bit) {
    return t & 1u << bit;
}

int main()
{
    std::chrono::seconds duration (60);
    std::srand(std::time(NULL));
    unsigned int i = std::rand();
    std::size_t bit = std::rand() % (sizeof(i) * CHAR_BIT);
    timer::default_alpha = 0.1l;
    std::cout
        << "func1 vs. func2: "
        << timer::compare(&func1, &func2, duration, i, bit)
        << std::endl;

    return 0;
}
