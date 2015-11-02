Timer
=====

A C++ header library for empirically estimating the time complexity of a
function, comparing 2 functions in terms, or just getting the raw timing for a
function.

Time Complexity Estimation
--------------------------

This calls the given function with increasing values of n, and then for each
time complexity that the timer can detect (currently log n, n, n log n, n^2,
n^3, n^4, 2^n, n!), it runs a statistical analysis to find the strength of the
association between the time complexity function and how long the given
function took in practice.

For example, try calling

    timer::approx_time_complexity([] (size_t n) {
        std::vector<size_t> v (n);
        for (size_t i = 0; i < n; ++i) {
            v[i] = n - i;
        }
        std::sort(v.begin(), v.end());
    }, 10000000, std::chrono::seconds(60));

and check out the return value, which should be equal to
`timer::timer_complextity::LINEARITHMIC` (AKA O(n log n)). This call should
last last pretty close to 60 seconds. The 10000000 means that the lambda
function won't be called with n > 10000000.

What's with timer.c and timer.h?
--------------------------------

These are a start at a pure C interface for this, with more features, but
they're not at all complete yet.
