#define BOOST
#include "timer.h"

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
    // auto func1 = [] () { };
    // auto func2 = [] () { };
    std::chrono::seconds duration (60);
    std::srand(std::time(NULL));
    unsigned int i = std::rand();
    std::size_t bit = std::rand() % (sizeof(i) * CHAR_BIT);
    timer::default_alpha = 0.1l;
    std::cout
        << "func2 vs. func3: "
        << timer::compare(&func2, &func3, duration, i, bit)
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
