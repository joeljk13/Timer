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
    unsigned int reps = 1024 * 1024;
    std::chrono::seconds dur (60);

    timer tm1 (func1, reps);
    assert(tm1.get_repetitions() == reps);
    std::cout << tm1.get_ratio() << std::endl;
    std::cout << "This should take 1 minute" << std::endl;
    timer tm2 (func2, dur);

    std::cout << tm1.get_ratio()
              << std::endl
              << tm2.get_ratio()
              << std::endl;

    return 0;
}
