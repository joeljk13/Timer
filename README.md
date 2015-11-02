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

Comparing Functions
-------------------

This times both functions and runs a statistical analysis (i.e. student's t
test) to determine which function is faster by a statistically significant
amount.

For example, if you want to compare sorting vectors vs. sorting deques, try
calling

    timer::compare([] () {
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

with `constexpr size_t n = 1000000` (or whatever number you want). It will run
for around 60 seconds, and then return 1 if the deque is faster, 0 if neither
is significantly faster than the other (NOT that they're statistically the
same), or -1 if the vector is faster (When I ran this, I got -1).

What's with timer.c and timer.h?
--------------------------------

These are a start at a pure C interface for this, with more features, but
they're not at all complete yet.
