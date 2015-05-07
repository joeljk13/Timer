#ifndef TIMER_HPP_
#define TIMER_HPP_ 1

#include <algorithm>
#include <chrono>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#include <cmath>
#include <cstddef>

class timer {

public:

    typedef std::chrono::nanoseconds     duration_type;
    typedef unsigned int                 repetition_type;
    typedef long long                    duration_per_repetition_type;
    typedef duration_per_repetition_type ratio_type;
    typedef std::chrono::steady_clock    clock_type;

    timer() :
        duration_   (0),
        reps_       (0)
    { }

    template <class Func, class... Args>
    timer(Func func, repetition_type reps, Args... args) :
        duration_   (0),
        reps_       (0)
    {
        this->measure(func, reps, std::forward<Args>(args)...);
    }

    template <class Func, class Reps, class Period, class... Args>
    timer(Func func, std::chrono::duration<Reps, Period> duration,
          Args... args) :
        duration_   (0),
        reps_       (0)
    {
        this->measure(func, duration, std::forward<Args>(args)...);
    }

    template <class Func, class... Args>
    void measure(Func func, repetition_type reps, Args... args)
    {
        clock_type::time_point start, end;

        reps_ += reps;

        start = clock_type::now();
        while (reps) {
            --reps;
            func(std::forward<Args>(args)...);
        }
        end = clock_type::now();

        duration_ += std::chrono::duration_cast<duration_type>(end - start);
    }

    template <class Func, class Reps, class Period, class... Args>
    void measure(Func func, std::chrono::duration<Reps, Period> duration,
                 Args... args)
    {
        clock_type::time_point start, end,
                               func_start (clock_type::now()),
                               func_end   (func_start + duration);

        duration_type           func_duration (0);
        constexpr duration_type zero = duration_type::zero();
        repetition_type         min_reps;

        for (min_reps = 1; func_duration != zero; min_reps *= 2) {
            repetition_type reps = min_reps;

            start = clock_type::now();
            while (reps) {
                --reps;
                func(std::forward<Args>(args)...);
            }
            end = clock_type::now();

            func_duration = std::chrono::duration_cast<duration_type>(end -
                start);
        }

        reps_     += min_reps;
        duration_ += func_duration;

        for (duration_type dur;
             (dur = std::chrono::duration_cast<duration_type>(func_end - end))
                > zero; ) {

            // Try to have this set of repetions spend half of the remaining
            // time alloted; while this will have more overhead, there's less
            // chance that it will accidentally go over the time limit
            repetition_type reps = static_cast<repetition_type>((dur * reps_) /
                (duration_ * 2));

            if (reps < min_reps) {
                return;
            }

            reps_ += reps;

            start = clock_type::now();
            while (reps) {
                --reps;
                func(std::forward<Args>(args)...);
            }
            end = clock_type::now();

            duration_ += std::chrono::duration_cast<duration_type>(end -
                start);
        }
    }

    duration_type get_duration() const
    {
        return duration_;
    }

    repetition_type get_repetitions() const
    {
        return reps_;
    }

    ratio_type get_duration_per_repetition() const
    {
        return duration_.count() / reps_;
    }

    ratio_type get_ratio() const
    {
        return get_duration_per_repetition();
    }

    typedef long double real;

    static constexpr real default_alpha = 0.1l;

    template <template <class, class...> class container>
    static int compare(container<timer> timers1, container<timer> timers2,
                       real alpha = default_alpha);

    static constexpr repetition_type MIN_REPS   = 8;
    static constexpr unsigned int    MAX_TIMERS = 16;

    template <class Func1, class Func2, class... Args>
    static int compare(Func1 func1, Func2 func2, repetition_type reps,
                       Args... args);

    template <class Func1, class Func2, class Reps, class Period,
              class... Args>
    static int compare(Func1 func1, Func2 func2,
                       std::chrono::duration<Reps, Period> duration,
                       Args... args);

    template <class T, std::size_t N>
    static real cc(T (&xs)[N], T (&ys)[N]);

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
    static time_complexity approx_time_complexity(Func func, std::size_t max_n,
        T duration_);

private:

    duration_type duration_;

    repetition_type reps_;

}; // class timer

template <template <class, class...> class container>
int timer::compare(container<timer> timers1, container<timer> timers2,
                   real alpha)
{
    std::vector<real> ratios1, ratios2;
    int size1 = 0,
        size2 = 0;

    real mean1 = 0.0l,
         mean2 = 0.0l,
         variance1 = 0.0l,
         variance2 = 0.0l,
         residuals = 0.0l;

    for (const timer &tm : timers1) {
        real ratio = static_cast<real>(tm.get_duration_per_repetition());
        mean1 += ratio;
        ++size1;
        ratios1.emplace_back(ratio);
    }

    for (const timer & tm : timers2) {
        real ratio = static_cast<real>(tm.get_duration_per_repetition());
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
    variance1 = (variance1 - residuals * residuals / size1) / (size1 - 1);
    residuals = 0.0l;

    for (real ratio : ratios2) {
        residuals += ratio;
        variance2 += ratio * ratio;
    }

    // Corrected 2-pass formula
    variance2 = (variance2 - residuals * residuals / size2) / (size2 - 1);
    real median1 = ratios1[size1 / 2];

    if (size1 % 2 == 0) {
        median1 = (median1 + ratios1[size1 / 2 + 1]) * 0.5l;
    }

    real median2 = ratios2[size2 / 2];
    if (size2 % 2 == 0) {
        median2 = (median2 + ratios2[size2 / 2 + 1]) * 0.5l;
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

    real t = (mean1 - mean2) / std::sqrt(tmp1 + tmp2);

    real df = (tmp1 + tmp2) * (tmp1 + tmp2) / (tmp1 * tmp1 / (size1 - 1) + tmp2
        * tmp2 / (size2 - 1));

#ifdef BOOST_VERSION

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
    real x = df / (t * t + df),
         a = df * 0.5l,
         b = 0.5l,
         tmp = (x <= 0.0l || x >= 1.0l ? 0.0l : std::exp(std::lgamma(a + b) + a
                 * std::log(x) - std::lgamma(a) + b * std::log(1.0l - x) -
                 std::lgamma(b)));

    bool sub_from_1 = false;

    if (x >= (a + 1.0l) / (a + b + 2.0l)) {
        sub_from_1 = true;
        std::swap(a, b);
        x = 1.0l - x;
    }

    // Use Lentz's method with continued fractions
    // I don't have good variable names, so these are just single-letter
    real c = 1.0l,
         d = 1.0l / (1.0l - (a + b) * x / (a + 1.0l)), e;

    // This doesn't represent the probability yet, but it will
    real prob = d;
    for (int n = 1; n <= 64; ++n) {
        e = n * (b - n) * x / ((a - 1.0l + n * 2) * (a + n * 2));
        d = 1.0l / (e * d + 1.0l);
        c = e / c + 1.0l;
        prob *= c * d;
        e = -(a + n) * (a + b + n) * x / ((a + n * 2) * (a + 1.0l + n * 2));
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
int timer::compare(Func1 func1, Func2 func2, timer::repetition_type reps,
                   Args... args)
{
    timer::repetition_type n = 1;
    std::vector<timer> timers1, timers2;

    while (reps / (n * 2) >= timer::MIN_REPS
        && n * 2 <= timer::MAX_TIMERS) {

        n *= 2;
    }

    reps /= n;
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

template <class Func1, class Func2, class Reps, class Period, class... Args>
int timer::compare(Func1 func1, Func2 func2,
                   std::chrono::duration<Reps, Period> duration,
                   Args... args)
{
    // Add 1 to allow for the time involved in the analysis
    duration_type dur (std::chrono::duration_cast<duration_type>(duration) /
        (2 * MAX_TIMERS + 1));

    std::vector<timer> timers1, timers2;

    timers1.reserve(MAX_TIMERS);
    timers2.reserve(MAX_TIMERS);
    repetition_type n = MAX_TIMERS / 2;

    while (n) {
        --n;
        timers1.emplace_back(func1, dur,
                std::forward<Args>(args)...);
        timers2.emplace_back(func2, dur,
                std::forward<Args>(args)...);
        timers2.emplace_back(func2, dur,
                std::forward<Args>(args)...);
        timers1.emplace_back(func1, dur,
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
    ratio_type times[N],
               dur_orig[N],
               dur[N];

    for (std::size_t i = 0, n = 1; i < N; ++i, n += max_n / N) {
        dur_orig[i] = n;
        timer tm (func, duration / N, n);
        times[i] = tm.get_ratio();
    }

    std::vector< std::pair<time_complexity, real> > crs;
    crs.reserve(N);

    // O(1) is easy to detect manually, plus it's difficult to detect it via
    // correlation, so just assume it's not that and move on

    // Logarithmic
    for (ratio_type i = 0; i < N; ++i) {
        dur[i] = std::ilogb(dur_orig[i]);
    }
    crs.emplace_back(std::make_pair(LOGARITHMIC, cc(times, dur)));

    // Linear
    for (ratio_type i = 0; i < N; ++i) {
        dur[i] = dur_orig[i];
    }
    crs.emplace_back(std::make_pair(LINEAR, cc(times, dur)));

    // Linearithmic
    for (ratio_type i = 0; i < N; ++i) {
        dur[i] = std::ilogb(dur_orig[i]) * dur_orig[i];
    }
    crs.emplace_back(std::make_pair(LINEARITHMIC, cc(times, dur)));

    // Quadratic
    for (ratio_type i = 0; i < N; ++i) {
        dur[i] = dur_orig[i] * dur_orig[i];
    }
    crs.emplace_back(std::make_pair(QUADRATIC, cc(times, dur)));

    // Cubic
    for (ratio_type i = 0; i < N; ++i) {
        ratio_type tmp = dur_orig[i];
        dur[i] = tmp * tmp * tmp;
    }
    crs.emplace_back(std::make_pair(CUBIC, cc(times, dur)));

    // Quartic
    for (ratio_type i = 0; i < N; ++i) {
        ratio_type tmp = dur_orig[i];
        dur[i] = (tmp * tmp) * (tmp * tmp);
    }
    crs.emplace_back(std::make_pair(QUADRATIC, cc(times, dur)));

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

            dur[i] = pow2(dur_orig[i]);
        }

        crs.emplace_back(std::make_pair(EXPONENTIAL, cc(times, dur)));
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

            dur[i] = factorial(dur_orig[i]);
        }

        crs.emplace_back(std::make_pair(FACTORIAL, cc(times, dur)));
    }

    auto iter = std::max_element(crs.begin(), crs.end(),
        [] (const std::pair<time_complexity, real> & elem1,
            const std::pair<time_complexity, real> & elem2) {

            return elem1.second < elem2.second;
        });

    return iter->first;
}

#endif
