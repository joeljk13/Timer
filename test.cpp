#include "timer.hpp"

#include <chrono>
#include <iostream>

#include <cassert>

static void func1()
{
    std::vector<int> v = {5, 3, 1, 3, 5, 3, 1, 3, 5, 3, 1, 3, 5, 3, 1, 3, 5};
    std::sort(v.begin(), v.end());
}

static void func2()
{
    std::vector<int> v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8};
    std::sort(v.begin(), v.end());
}

int main(void)
{
    unsigned int reps = 1024 * 1024;
    std::chrono::seconds dur (60);

    std::cout << "This should take 1 minute" << std::endl;

    timer tm1 (func1, reps);
    assert(tm1.get_repetitions() == reps);
    timer tm2 (func2, dur);

    std::cout << tm1.get_ratio()
              << std::endl
              << tm2.get_ratio()
              << std::endl;

    return 0;
}
