#include "timer.hpp"

#include <chrono>
#include <deque>
#include <iostream>
#include <vector>

int main()
{
    // Time complexity example
    timer::time_complexity tc = timer::approx_time_complexity([] (size_t n) {
        std::vector<size_t> v (n);
        for (size_t i = 0; i < n; ++i) {
            v[i] = n - i;
        }
        std::sort(v.begin(), v.end());
    }, 10000000, std::chrono::seconds(60));

    std::cout << "This should be " << (int)timer::time_complexity::LINEARITHMIC << ":\n";
    std::cout << (int)tc << std::endl;

    // Function comparison example - vector vs. deque
    constexpr size_t n = 1000000;

    std::cout << timer::compare([] () {
        std::deque<int> d (n);
        for (size_t i = 0; i < n; ++i) {
           d[i] = n - i;
        }
        std::sort(d.begin(), d.end());
    }, [] () {
        std::vector<int> v (n);
        for (size_t i = 0; i < n; ++i) {
            v[i] = n - i;
        }
        std::sort(v.begin(), v.end());
    }, std::chrono::seconds(60));

    return 0;
}
