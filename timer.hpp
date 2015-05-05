#ifndef h_
#define h_

// If HAVE_BOOST or _HAVE_BOOST or BOOST or _BOOST is defined, this uses the
// boost library for the statistical analysis; otherwise, it does some built in
// math to get a good enough answer. Boost is preferred, if possible.

#include <algorithm>
#include <chrono>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#if !defined(BOOST) && \
    (defined(HAVE_BOOST)   || \
    defined(_HAVE_BOOST)   || \
    defined(BOOST)         || \
    defined(_BOOST))
#define BOOST
#endif

#include <cmath>
#include <cstddef>

#if defined(BOOST)
#include <boost/math/distributions/students_t.hpp>
#endif

namespace _timer = std::chrono;

class timer {

public:

    typedef _timer::nanoseconds duration_type;
    typedef unsigned int repetition_type;
    typedef long long duration_per_repetition_type;
    typedef duration_per_repetition_type ratio_type;
    typedef _timer::steady_clock clock_type;

    timer() :
        duration_ (0),
        reps_ (0)
    { }

    template <class Func, class ... Args>
    timer(Func func, repetition_type reps,
            Args... args) :
        duration_ (0),
        reps_ (0)
    {
        this->measure(func, reps, std::forward<Args>(
                args)...);
    }

    template <class Func, class Reps, class Period,
            class... Args>
    timer(Func func, _timer::duration<Reps, Period>
            duration, Args... args) :
        duration_ (0),
        reps_ (0)
    {
        this->measure(func, duration, std::forward<Args>(
                args)...);
    }

    template <class Func, class... Args>
    timer & measure(Func func, repetition_type reps,
            Args... args)
    {
        reps_ += reps;
        clock_type::time_point start, end;
        start = clock_type::now();
        while (reps) {
            --reps;
            func(std::forward<Args>(args)...);
        }
        end = clock_type::now();
        duration_ += _timer::duration_cast<duration_type>(end -
                start);
        return *this;
    }

    template <class Func, class Reps, class Period,
            class... Args>
    timer & measure(Func func, _timer::duration<Reps,
            Period> duration, Args... args)
    {
        clock_type::time_point func_start (clock_type::now());
        clock_type::time_point start, end, func_end
                (func_start + duration);
        constexpr duration_type zero = duration_type::zero();
        duration_type func_duration (0);
        repetition_type r = 1;
        while (true) {
            repetition_type reps = r;
            start = clock_type::now();
            while (reps) {
                --reps;
                func(std::forward<Args>(args)...);
            }
            end = clock_type::now();
            func_duration = _timer::duration_cast<duration_type>(
                    end - start);
            if (func_duration != zero) {
                break;
            }
            r *= 2;
        }
        reps_ += r;
        duration_ += func_duration;
        duration_type ns;
        while ((ns = _timer::duration_cast<duration_type>(
                func_end - end)) > zero) {
            // Try to have this set of repetions spend half of the remaining
            // time alloted; while this will have more overhead, there's less
            // chance that it will accidentally go over the time limit
            repetition_type reps = static_cast<repetition_type>(
                    (ns * reps_) / (duration_ * 2));
            if (reps < r) {
                return *this;
            }
            reps_ += reps;
            start = clock_type::now();
            while (reps) {
                --reps;
                func(std::forward<Args>(args)...);
            }
            end = clock_type::now();
            duration_ += _timer::duration_cast<duration_type>(
                    end - start);
        }
        return *this;
    }

    duration_type get_duration() const
    { return duration_; }

    repetition_type get_repetitions() const
    { return reps_; }

    ratio_type getduration_per_repetition() const
    { return duration_.count() / reps_; }

    ratio_type get_ratio() const
    { return getduration_per_repetition(); }

    typedef long double real;

    static real default_alpha;

    template <template <class, class...> class C>
    static int compare(C<timer> timers1, C<timer> timers2,
            real alpha = default_alpha);

    static constexpr repetition_type MIN_REPS = 8;
    static constexpr unsigned int MAX_TIMERS = 16;

    template <class Func1, class Func2, class... Args>
    static int compare(Func1 func1, Func2 func2,
            repetition_type reps, Args... args);

    template <class Func1, class Func2, class Reps,
            class Period, class... Args>
    static int compare(Func1 func1, Func2 func2,
            _timer::duration<Reps, Period> duration,
            Args... args);

    template <class T, std::size_t N>
    static real cc(T (&xs)[N],
                                 T (&ys)[N]);

    enum time_complexity {
        CONSTANT,
        LOGARITHMIC,
        LINEAR,
        LINEARITHMIC,
        QUADRATIC,
        CUBIC,
        QUARTIC,
        EXPONENTIAL,
        FACTORIAL
    };

    template <class Func, class T, std::size_t N = 16>
    static time_complexity approx_time_complexity(Func func,
            std::size_t max_n, T duration_);

private:

    duration_type duration_;

    repetition_type reps_;

}; // class timer

timer::real timer::default_alpha = 0.01l;

template <template <class, class...> class C>
int timer::compare(C<timer> timers1, C<timer> timers2,
        real alpha)
{
    std::vector<real> ratios1, ratios2;
    int size1 = 0, size2 = 0;
    real mean1 = 0.0l, mean2 = 0.0l,
            variance1 = 0.0l, variance2 = 0.0l,
            residuals = 0.0l;
    for (const timer & tm : timers1) {
        real ratio = static_cast<real>(
                tm.getduration_per_repetition());
        mean1 += ratio;
        ++size1;
        ratios1.emplace_back(ratio);
    }
    for (const timer & tm : timers2) {
        real ratio = static_cast<real>(
                tm.getduration_per_repetition());
        mean2 += ratio;
        ++size2;
        ratios2.emplace_back(ratio);
    }
    mean1 /= size1;
    mean2 /= size2;
    std::sort(ratios1.begin(), ratios1.end());
    std::sort(ratios2.begin(), ratios2.end());
    for (real ratio : ratios1) {
        residuals += ratio;
        variance1 += ratio * ratio;
    }
    // Corrected 2-pass formula
    variance1 = (variance1 - residuals * residuals
            / size1) / (size1 - 1);
    residuals = 0.0l;
    for (real ratio : ratios2) {
        residuals += ratio;
        variance2 += ratio * ratio;
    }
    // Corrected 2-pass formula
    variance2 = (variance2 - residuals * residuals
            / size2) / (size2 - 1);
    real median1 = ratios1[size1 / 2];
    if (size1 % 2 == 0) {
        median1 = (median1 + ratios1[size1 / 2 +
                1]) * 0.5l;
    }
    real median2 = ratios2[size2 / 2];
    if (size2 % 2 == 0) {
        median2 = (median2 + ratios2[size2 / 2 +
                1]) * 0.5l;
    }
    // Test that timers1 < timers2
    bool swapped = false;
    if (mean2 < mean1 && median2 < median1) {
        swapped = true;
        // Only swap what will be used later on
        std::swap(size1, size2);
        std::swap(mean1, mean2);
        std::swap(median1, median2);
        std::swap(variance1, variance2);
    }
    if (mean1 >= mean2 || median1 >= median2) {
        // Neither set won in both mean and median
        return 0;
    }
    // t-test variables
    real tmp1 = variance1 / size1,
                tmp2 = variance2 / size2;
    real t = (mean1 - mean2) /
            std::sqrt(tmp1 + tmp2);
    real df = (tmp1 + tmp2) *
                            (tmp1 + tmp2) /
                            (tmp1 * tmp1 / (size1 - 1) +
                            tmp2 * tmp2 / (size2 - 1));
#if defined(BOOST)
    boost::math::students_t dist (static_cast<double>(df));
    real prob = boost::math::cdf(dist, t);
    if (prob < alpha) {
        int ret = (mean1 < mean2 ? 1 : -1);
        if (swapped) {
            ret = -ret;
        }
        return ret;
    }
    return 0;
#else
    // Incomplete beta function
    real x = df / (t * t + df);
    real a = df * 0.5l;
    real b = 0.5l;
    real tmp = (x <= 0.0l || x >= 1.0l ? 0.0l :
            std::exp(std::lgamma(a + b) +
            a * std::log(x) - std::lgamma(a) +
            b * std::log(1.0l - x) - std::lgamma(b)));
    bool sub_from_1 = false;
    if (x >= (a + 1.0l) / (a + b + 2.0l)) {
        sub_from_1 = true;
        std::swap(a, b);
        x = 1.0l - x;
    }
    // Use Lentz's method with continued fractions
    // I don't have good variable names, so these are just single-letter
    real c = 1.0l,
            d = 1.0l / (1.0l - (a + b) * x /
            (a + 1.0l)), e;
    // This doesn't represent the probability yet, but it will
    real prob = d;
    for (int n = 1; n <= 64; ++n) {
        e = n * (b - n) * x /
                ((a - 1.0l + n * 2) * (a + n * 2));
        d = 1.0l / (e * d + 1.0l);
        c = e / c + 1.0l;
        prob *= c * d;
        e = -(a + n) * (a + b + n) *
                x / ((a + n * 2) *
                (a + 1.0l + n * 2));
        d = 1.0l / (e * d + 1.0l);
        c = e / c + 1.0l;
        prob *= c * d;
        if (c * d - 1.0l < 1.0e-12l) {
            break;
        }
    }
    prob *= tmp / a;
    if (sub_from_1) {
        prob = 1.0l - prob;
    }
    if (prob < 0.0l || prob > 1.0l) {
        return 0;
    }
    if (prob < alpha / 2) {
        int ret = (mean1 < mean2 ? 1 : -1);
        if (swapped) {
            ret = -ret;
        }
        return ret;
    }
    return 0;
#endif
}

template <class Func1, class Func2, class... Args>
int timer::compare(Func1 func1, Func2 func2,
        timer::repetition_type reps, Args... args)
{
    timer::repetition_type n = 1;
    while (reps / (n * 2) >= timer::MIN_REPS &&
            n * 2 <= timer::MAX_TIMERS) {
        n *= 2;
    }
    reps /= n;
    std::vector<timer> timers1, timers2;
    timers1.reserve(n);
    timers2.reserve(n);
    n /= 2;
    while (n) {
        --n;
        timers1.emplace_back(func1, reps,
                std::forward<Args>(args)...);
        timers2.emplace_back(func2, reps,
                std::forward<Args>(args)...);
        timers2.emplace_back(func2, reps,
                std::forward<Args>(args)...);
        timers1.emplace_back(func1, reps,
                std::forward<Args>(args)...);
    }
    return compare(timers1, timers2);
}

template <class Func1, class Func2, class Reps,
        class Period, class... Args>
int timer::compare(Func1 func1, Func2 func2,
        _timer::duration<Reps, Period> duration,
        Args... args)
{
    // Add 1 to allow for the time involved in the analysis
    duration_type ns
            (_timer::duration_cast<duration_type>(duration) /
            (2 * MAX_TIMERS + 1));
    std::vector<timer> timers1, timers2;
    timers1.reserve(MAX_TIMERS);
    timers2.reserve(MAX_TIMERS);
    repetition_type n = MAX_TIMERS / 2;
    while (n) {
        --n;
        timers1.emplace_back(func1, ns,
                std::forward<Args>(args)...);
        timers2.emplace_back(func2, ns,
                std::forward<Args>(args)...);
        timers2.emplace_back(func2, ns,
                std::forward<Args>(args)...);
        timers1.emplace_back(func1, ns,
                std::forward<Args>(args)...);
    }
    return compare(timers1, timers2);
}

// Calculates the correlation coefficient
template <class T, std::size_t N>
timer::real timer::cc(T (&xs)[N],
                                    T (&ys)[N])
{
    real x[N];
    real y[N];
    real xmean = 0.0l, ymean = 0.0l;
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = static_cast<real>(xs[i]);
        xmean += x[i];
        y[i] = static_cast<real>(ys[i]);
        ymean += y[i];
    }
    xmean /= N;
    ymean /= N;
    real xres, yres, sx = 0.0l, sy = 0.0l,
            sxy = 0.0l;
    for (std::size_t i = 0; i < N; ++i) {
        xres = x[i] - xmean;
        sx += xres * xres;
        yres = y[i] - ymean;
        sy += yres * yres;
        sxy += xres * yres;
    }
    return sxy / std::sqrt(sx * sy);
}

template <class Func, class T, std::size_t N>
timer::time_complexity timer::approx_time_complexity(Func func,
        std::size_t max_n, T duration)
{
    ratio_type times[N];
    ratio_type ns_orig[N];
    for (std::size_t i = 0, n = 1;
            i < N;
            ++i, n += max_n / N) {
        ns_orig[i] = n;
        timer tm (func, duration / N, n);
        times[i] = tm.get_ratio();
    }
    std::vector< std::pair<time_complexity, real> > crs;
    crs.reserve(N);
    // O(1) is easy to detect manually, plus it's difficult to detect it via
    // correlation, so just assume it's not that and move on
    ratio_type ns[N];
    // Logarithmic
    for (ratio_type i = 0; i < N; ++i) {
        ns[i] = std::ilogb(ns_orig[i]);
    }
    crs.emplace_back(std::make_pair(LOGARITHMIC,
            cc(times, ns)));
    // Linear
    for (ratio_type i = 0; i < N; ++i) {
        ns[i] = ns_orig[i];
    }
    crs.emplace_back(std::make_pair(LINEAR,
            cc(times, ns)));
    // Linearithmic
    for (ratio_type i = 0; i < N; ++i) {
        ns[i] = std::ilogb(ns_orig[i]) *
                ns_orig[i];
    }
    crs.emplace_back(std::make_pair(LINEARITHMIC,
            cc(times, ns)));
    // Quadratic
    for (ratio_type i = 0; i < N; ++i) {
        ns[i] = ns_orig[i] *
                ns_orig[i];
    }
    crs.emplace_back(std::make_pair(QUADRATIC,
            cc(times, ns)));
    // Cubic
    for (ratio_type i = 0; i < N; ++i) {
        ratio_type tmp = ns_orig[i];
        ns[i] = tmp * tmp * tmp;
    }
    crs.emplace_back(std::make_pair(CUBIC,
            cc(times, ns)));
    // Quartic
    for (ratio_type i = 0; i < N; ++i) {
        ratio_type tmp = ns_orig[i];
        ns[i] = (tmp * tmp) *
                (tmp * tmp);
    }
    crs.emplace_back(std::make_pair(QUADRATIC,
            cc(times, ns)));
    // Only calclate these next ones if they won't overflow
    // Exponential
    if (max_n < std::ilogb(std::numeric_limits<ratio_type>::max())) {
        for (ratio_type i = 0; i < N; ++i) {
            auto pow2 = [] (ratio_type e) {
                ratio_type result = 1;
                ratio_type b = 2;
                for ( ; e > 0; e /= 2) {
                    if (e % 2 == 1) {
                        result *= b;
                    }
                    b *= b;
                }
                return result;
            };
            ns[i] = pow2(ns_orig[i]);
        }
        crs.emplace_back(std::make_pair(EXPONENTIAL,
                cc(times, ns)));
    }
    // Factorial
    if (max_n <= 20) {
        for (ratio_type i = 0; i < N; ++i) {
            auto factorial = [] (ratio_type x) -> ratio_type {
                ratio_type f = 1;
                for (ratio_type j = 2; j <= x; ++j) {
                    f *= j;
                }
                return f;
            };
            ns[i] = factorial(ns_orig[i]);
        }
        crs.emplace_back(std::make_pair(FACTORIAL,
                cc(times, ns)));
    }
    auto iter = std::max_element(crs.begin(), crs.end(),
            [] (const std::pair<time_complexity, real> & elem1,
                const std::pair<time_complexity, real> & elem2)
            {
                return elem1.second < elem2.second;
            });
    return iter->first;
}

#endif // h_
