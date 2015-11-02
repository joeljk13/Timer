#include "timer.hpp"

#include <chrono>
#include <iostream>

#include <cassert>

static void func1()
{
    std::vector<int> v;

    for (int i = 0; i < 100; ++i) {
        v.push_back(i);
    }

    std::sort(v.begin(), v.end());
}

static void func2()
{
    std::vector<int> v;

    for (int i = 0; i < 128; ++i) {
        v.push_back(i);
    }

    std::sort(v.begin(), v.end());
}

int main(void)
{
    timer::time_complexity tc = timer::approx_time_complexity([] (size_t n) {
        std::vector<size_t> v (n);
        for (size_t i = 0; i < n; ++i) {
            v[i] = n - i;
        }
        std::sort(v.begin(), v.end());
    }, 10000000, std::chrono::seconds(60));

    std::cout << (int)tc << std::endl;
    std::cout << (int)timer::time_complexity::LINEARITHMIC << std::endl;

    return 0;
}
