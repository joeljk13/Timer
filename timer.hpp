#ifndef _timer_h_
#define _timer_h_

// All non-public names begin with _timer to avoid naming conflicts

// If HAVE_BOOST or _HAVE_BOOST or BOOST or _BOOST is defined, this uses the
// boost library for the statistical analysis; otherwise, it does some built in
// math to get a good enough answer. Boost is preferred, if possible.

#include <algorithm>
#include <chrono>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#if !defined(_timer_BOOST) && \
    (defined(HAVE_BOOST)   || \
    defined(_HAVE_BOOST)   || \
    defined(BOOST)         || \
    defined(_BOOST))
#define _timer_BOOST
#endif

#include <cmath>
#include <cstddef>

#if defined(_timer_BOOST)
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
        _timer_duration_ (0),
        _timer_reps_ (0)
    { }

    template <class _timer_Func, class ... _timer_Args>
    timer(_timer_Func _timer_func, repetition_type _timer_reps,
            _timer_Args... _timer_args) :
        _timer_duration_ (0),
        _timer_reps_ (0)
    {
        this->measure(_timer_func, _timer_reps, std::forward<_timer_Args>(
                _timer_args)...);
    }

    template <class _timer_Func, class _timer_Reps, class _timer_Period,
            class... _timer_Args>
    timer(_timer_Func _timer_func, _timer::duration<_timer_Reps, _timer_Period>
            _timer_duration, _timer_Args... _timer_args) :
        _timer_duration_ (0),
        _timer_reps_ (0)
    {
        this->measure(_timer_func, _timer_duration, std::forward<_timer_Args>(
                _timer_args)...);
    }

    template <class _timer_Func, class... _timer_Args>
    timer & measure(_timer_Func _timer_func, repetition_type _timer_reps,
            _timer_Args... _timer_args)
    {
        _timer_reps_ += _timer_reps;
        clock_type::time_point _timer_start, _timer_end;
        _timer_start = clock_type::now();
        while (_timer_reps) {
            --_timer_reps;
            _timer_func(std::forward<_timer_Args>(_timer_args)...);
        }
        _timer_end = clock_type::now();
        _timer_duration_ += _timer::duration_cast<duration_type>(_timer_end -
                _timer_start);
        return *this;
    }

    template <class _timer_Func, class _timer_Reps, class _timer_Period,
            class... _timer_Args>
    timer & measure(_timer_Func _timer_func, _timer::duration<_timer_Reps,
            _timer_Period> _timer_duration, _timer_Args... _timer_args)
    {
        clock_type::time_point _timer_func_start (clock_type::now());
        clock_type::time_point _timer_start, _timer_end, _timer_func_end
                (_timer_func_start + _timer_duration);
        constexpr duration_type _timer_zero = duration_type::zero();
        duration_type _timer_func_duration (0);
        repetition_type _timer_r = 1;
        while (true) {
            repetition_type _timer_reps = _timer_r;
            _timer_start = clock_type::now();
            while (_timer_reps) {
                --_timer_reps;
                _timer_func(std::forward<_timer_Args>(_timer_args)...);
            }
            _timer_end = clock_type::now();
            _timer_func_duration = _timer::duration_cast<duration_type>(
                    _timer_end - _timer_start);
            if (_timer_func_duration != _timer_zero) {
                break;
            }
            _timer_r *= 2;
        }
        _timer_reps_ += _timer_r;
        _timer_duration_ += _timer_func_duration;
        duration_type _timer_ns;
        while ((_timer_ns = _timer::duration_cast<duration_type>(
                _timer_func_end - _timer_end)) > _timer_zero) {
            // Try to have this set of repetions spend half of the remaining
            // time alloted; while this will have more overhead, there's less
            // chance that it will accidentally go over the time limit
            repetition_type _timer_reps = static_cast<repetition_type>(
                    (_timer_ns * _timer_reps_) / (_timer_duration_ * 2));
            if (_timer_reps < _timer_r) {
                return *this;
            }
            _timer_reps_ += _timer_reps;
            _timer_start = clock_type::now();
            while (_timer_reps) {
                --_timer_reps;
                _timer_func(std::forward<_timer_Args>(_timer_args)...);
            }
            _timer_end = clock_type::now();
            _timer_duration_ += _timer::duration_cast<duration_type>(
                    _timer_end - _timer_start);
        }
        return *this;
    }

    duration_type get_duration() const
    { return _timer_duration_; }

    repetition_type get_repetitions() const
    { return _timer_reps_; }

    ratio_type get_timer_duration_per_repetition() const
    { return _timer_duration_.count() / _timer_reps_; }

    ratio_type get_ratio() const
    { return get_timer_duration_per_repetition(); }

    typedef long double _timer_real;

    static _timer_real default_alpha;

    template <template <class, class...> class _timer_C>
    static int compare(_timer_C<timer> _timer_timers1, _timer_C<timer> _timer_timers2,
            _timer_real _timer_alpha = default_alpha);

    static constexpr repetition_type MIN_REPS = 8;
    static constexpr unsigned int MAX_TIMERS = 16;

    template <class _timer_Func1, class _timer_Func2, class... _timer_Args>
    static int compare(_timer_Func1 _timer_func1, _timer_Func2 _timer_func2,
            repetition_type _timer_reps, _timer_Args... _timer_args);

    template <class _timer_Func1, class _timer_Func2, class _timer_Reps,
            class _timer_Period, class... _timer_Args>
    static int compare(_timer_Func1 _timer_func1, _timer_Func2 _timer_func2,
            _timer::duration<_timer_Reps, _timer_Period> _timer_duration,
            _timer_Args... _timer_args);

    template <class _timer_T, std::size_t _timer_N>
    static _timer_real _timer_cc(_timer_T (&_timer_xs)[_timer_N],
                                 _timer_T (&_timer_ys)[_timer_N]);

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

    template <class _timer_Func, class _timer_T, std::size_t _timer_N = 16>
    static time_complexity approx_time_complexity(_timer_Func _timer_func,
            std::size_t _timer_max_n, _timer_T _timer_duration_);

private:

    duration_type _timer_duration_;

    repetition_type _timer_reps_;

}; // class timer

timer::_timer_real timer::default_alpha = 0.01l;

template <template <class, class...> class _timer_C>
int timer::compare(_timer_C<timer> _timer_timers1, _timer_C<timer> _timer_timers2,
        _timer_real _timer_alpha)
{
    std::vector<_timer_real> _timer_ratios1, _timer_ratios2;
    int _timer_size1 = 0, _timer_size2 = 0;
    _timer_real _timer_mean1 = 0.0l, _timer_mean2 = 0.0l,
            _timer_variance1 = 0.0l, _timer_variance2 = 0.0l,
            _timer_residuals = 0.0l;
    for (const timer & _timer_tm : _timer_timers1) {
        _timer_real _timer_ratio = static_cast<_timer_real>(
                _timer_tm.get_timer_duration_per_repetition());
        _timer_mean1 += _timer_ratio;
        ++_timer_size1;
        _timer_ratios1.emplace_back(_timer_ratio);
    }
    for (const timer & _timer_tm : _timer_timers2) {
        _timer_real _timer_ratio = static_cast<_timer_real>(
                _timer_tm.get_timer_duration_per_repetition());
        _timer_mean2 += _timer_ratio;
        ++_timer_size2;
        _timer_ratios2.emplace_back(_timer_ratio);
    }
    _timer_mean1 /= _timer_size1;
    _timer_mean2 /= _timer_size2;
    std::sort(_timer_ratios1.begin(), _timer_ratios1.end());
    std::sort(_timer_ratios2.begin(), _timer_ratios2.end());
    for (_timer_real _timer_ratio : _timer_ratios1) {
        _timer_residuals += _timer_ratio;
        _timer_variance1 += _timer_ratio * _timer_ratio;
    }
    // Corrected 2-pass formula
    _timer_variance1 = (_timer_variance1 - _timer_residuals * _timer_residuals
            / _timer_size1) / (_timer_size1 - 1);
    _timer_residuals = 0.0l;
    for (_timer_real _timer_ratio : _timer_ratios2) {
        _timer_residuals += _timer_ratio;
        _timer_variance2 += _timer_ratio * _timer_ratio;
    }
    // Corrected 2-pass formula
    _timer_variance2 = (_timer_variance2 - _timer_residuals * _timer_residuals
            / _timer_size2) / (_timer_size2 - 1);
    _timer_real _timer_median1 = _timer_ratios1[_timer_size1 / 2];
    if (_timer_size1 % 2 == 0) {
        _timer_median1 = (_timer_median1 + _timer_ratios1[_timer_size1 / 2 +
                1]) * 0.5l;
    }
    _timer_real _timer_median2 = _timer_ratios2[_timer_size2 / 2];
    if (_timer_size2 % 2 == 0) {
        _timer_median2 = (_timer_median2 + _timer_ratios2[_timer_size2 / 2 +
                1]) * 0.5l;
    }
    // Test that timers1 < timers2
    bool _timer_swapped = false;
    if (_timer_mean2 < _timer_mean1 && _timer_median2 < _timer_median1) {
        _timer_swapped = true;
        // Only swap what will be used later on
        std::swap(_timer_size1, _timer_size2);
        std::swap(_timer_mean1, _timer_mean2);
        std::swap(_timer_median1, _timer_median2);
        std::swap(_timer_variance1, _timer_variance2);
    }
    if (_timer_mean1 >= _timer_mean2 || _timer_median1 >= _timer_median2) {
        // Neither set won in both mean and median
        return 0;
    }
    // t-test variables
    _timer_real _timer_tmp1 = _timer_variance1 / _timer_size1,
                _timer_tmp2 = _timer_variance2 / _timer_size2;
    _timer_real _timer_t = (_timer_mean1 - _timer_mean2) /
            std::sqrt(_timer_tmp1 + _timer_tmp2);
    _timer_real _timer_df = (_timer_tmp1 + _timer_tmp2) *
                            (_timer_tmp1 + _timer_tmp2) /
                            (_timer_tmp1 * _timer_tmp1 / (_timer_size1 - 1) +
                            _timer_tmp2 * _timer_tmp2 / (_timer_size2 - 1));
#if defined(_timer_BOOST)
    boost::math::students_t _timer_dist (static_cast<double>(_timer_df));
    _timer_real _timer_prob = boost::math::cdf(_timer_dist, _timer_t);
    if (_timer_prob < _timer_alpha) {
        int _timer_ret = (_timer_mean1 < _timer_mean2 ? 1 : -1);
        if (_timer_swapped) {
            _timer_ret = -_timer_ret;
        }
        return _timer_ret;
    }
    return 0;
#else
    // Incomplete beta function
    _timer_real _timer_x = _timer_df / (_timer_t * _timer_t + _timer_df);
    _timer_real _timer_a = _timer_df * 0.5l;
    _timer_real _timer_b = 0.5l;
    _timer_real _timer_tmp = (_timer_x <= 0.0l || _timer_x >= 1.0l ? 0.0l :
            std::exp(std::lgamma(_timer_a + _timer_b) +
            _timer_a * std::log(_timer_x) - std::lgamma(_timer_a) +
            _timer_b * std::log(1.0l - _timer_x) - std::lgamma(_timer_b)));
    bool _timer_sub_from_1 = false;
    if (_timer_x >= (_timer_a + 1.0l) / (_timer_a + _timer_b + 2.0l)) {
        _timer_sub_from_1 = true;
        std::swap(_timer_a, _timer_b);
        _timer_x = 1.0l - _timer_x;
    }
    // Use Lentz's method with continued fractions
    // I don't have good variable names, so these are just single-letter
    _timer_real _timer_c = 1.0l,
            _timer_d = 1.0l / (1.0l - (_timer_a + _timer_b) * _timer_x /
            (_timer_a + 1.0l)), _timer_e;
    // This doesn't represent the probability yet, but it will
    _timer_real _timer_prob = _timer_d;
    for (int _timer_n = 1; _timer_n <= 64; ++_timer_n) {
        _timer_e = _timer_n * (_timer_b - _timer_n) * _timer_x /
                ((_timer_a - 1.0l + _timer_n * 2) * (_timer_a + _timer_n * 2));
        _timer_d = 1.0l / (_timer_e * _timer_d + 1.0l);
        _timer_c = _timer_e / _timer_c + 1.0l;
        _timer_prob *= _timer_c * _timer_d;
        _timer_e = -(_timer_a + _timer_n) * (_timer_a + _timer_b + _timer_n) *
                _timer_x / ((_timer_a + _timer_n * 2) *
                (_timer_a + 1.0l + _timer_n * 2));
        _timer_d = 1.0l / (_timer_e * _timer_d + 1.0l);
        _timer_c = _timer_e / _timer_c + 1.0l;
        _timer_prob *= _timer_c * _timer_d;
        if (_timer_c * _timer_d - 1.0l < 1.0e-12l) {
            break;
        }
    }
    _timer_prob *= _timer_tmp / _timer_a;
    if (_timer_sub_from_1) {
        _timer_prob = 1.0l - _timer_prob;
    }
    if (_timer_prob < 0.0l || _timer_prob > 1.0l) {
        return 0;
    }
    if (_timer_prob < _timer_alpha / 2) {
        int _timer_ret = (_timer_mean1 < _timer_mean2 ? 1 : -1);
        if (_timer_swapped) {
            _timer_ret = -_timer_ret;
        }
        return _timer_ret;
    }
    return 0;
#endif
}

template <class _timer_Func1, class _timer_Func2, class... _timer_Args>
int timer::compare(_timer_Func1 _timer_func1, _timer_Func2 _timer_func2,
        timer::repetition_type _timer_reps, _timer_Args... _timer_args)
{
    timer::repetition_type _timer_n = 1;
    while (_timer_reps / (_timer_n * 2) >= timer::MIN_REPS &&
            _timer_n * 2 <= timer::MAX_TIMERS) {
        _timer_n *= 2;
    }
    _timer_reps /= _timer_n;
    std::vector<timer> _timer_timers1, _timer_timers2;
    _timer_timers1.reserve(_timer_n);
    _timer_timers2.reserve(_timer_n);
    _timer_n /= 2;
    while (_timer_n) {
        --_timer_n;
        _timer_timers1.emplace_back(_timer_func1, _timer_reps,
                std::forward<_timer_Args>(_timer_args)...);
        _timer_timers2.emplace_back(_timer_func2, _timer_reps,
                std::forward<_timer_Args>(_timer_args)...);
        _timer_timers2.emplace_back(_timer_func2, _timer_reps,
                std::forward<_timer_Args>(_timer_args)...);
        _timer_timers1.emplace_back(_timer_func1, _timer_reps,
                std::forward<_timer_Args>(_timer_args)...);
    }
    return compare(_timer_timers1, _timer_timers2);
}

template <class _timer_Func1, class _timer_Func2, class _timer_Reps,
        class _timer_Period, class... _timer_Args>
int timer::compare(_timer_Func1 _timer_func1, _timer_Func2 _timer_func2,
        _timer::duration<_timer_Reps, _timer_Period> _timer_duration,
        _timer_Args... _timer_args)
{
    // Add 1 to allow for the time involved in the analysis
    duration_type _timer_ns
            (_timer::duration_cast<duration_type>(_timer_duration) /
            (2 * MAX_TIMERS + 1));
    std::vector<timer> _timer_timers1, _timer_timers2;
    _timer_timers1.reserve(MAX_TIMERS);
    _timer_timers2.reserve(MAX_TIMERS);
    repetition_type _timer_n = MAX_TIMERS / 2;
    while (_timer_n) {
        --_timer_n;
        _timer_timers1.emplace_back(_timer_func1, _timer_ns,
                std::forward<_timer_Args>(_timer_args)...);
        _timer_timers2.emplace_back(_timer_func2, _timer_ns,
                std::forward<_timer_Args>(_timer_args)...);
        _timer_timers2.emplace_back(_timer_func2, _timer_ns,
                std::forward<_timer_Args>(_timer_args)...);
        _timer_timers1.emplace_back(_timer_func1, _timer_ns,
                std::forward<_timer_Args>(_timer_args)...);
    }
    return compare(_timer_timers1, _timer_timers2);
}

// Calculates the correlation coefficient
template <class _timer_T, std::size_t _timer_N>
timer::_timer_real timer::_timer_cc(_timer_T (&_timer_xs)[_timer_N],
                                    _timer_T (&_timer_ys)[_timer_N])
{
    _timer_real _timer_x[_timer_N];
    _timer_real _timer_y[_timer_N];
    _timer_real _timer_xmean = 0.0l, _timer_ymean = 0.0l;
    for (std::size_t _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
        _timer_x[_timer_i] = static_cast<_timer_real>(_timer_xs[_timer_i]);
        _timer_xmean += _timer_x[_timer_i];
        _timer_y[_timer_i] = static_cast<_timer_real>(_timer_ys[_timer_i]);
        _timer_ymean += _timer_y[_timer_i];
    }
    _timer_xmean /= _timer_N;
    _timer_ymean /= _timer_N;
    _timer_real _timer_xres, _timer_yres, _timer_sx = 0.0l, _timer_sy = 0.0l,
            _timer_sxy = 0.0l;
    for (std::size_t _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
        _timer_xres = _timer_x[_timer_i] - _timer_xmean;
        _timer_sx += _timer_xres * _timer_xres;
        _timer_yres = _timer_y[_timer_i] - _timer_ymean;
        _timer_sy += _timer_yres * _timer_yres;
        _timer_sxy += _timer_xres * _timer_yres;
    }
    return _timer_sxy / std::sqrt(_timer_sx * _timer_sy);
}

template <class _timer_Func, class _timer_T, std::size_t _timer_N>
timer::time_complexity timer::approx_time_complexity(_timer_Func _timer_func,
        std::size_t _timer_max_n, _timer_T _timer_duration)
{
    ratio_type _timer_times[_timer_N];
    ratio_type _timer_ns_orig[_timer_N];
    for (std::size_t _timer_i = 0, _timer_n = 1;
            _timer_i < _timer_N;
            ++_timer_i, _timer_n += _timer_max_n / _timer_N) {
        _timer_ns_orig[_timer_i] = _timer_n;
        timer _timer_tm (_timer_func, _timer_duration / _timer_N, _timer_n);
        _timer_times[_timer_i] = _timer_tm.get_ratio();
    }
    std::vector< std::pair<time_complexity, _timer_real> > _timer_crs;
    _timer_crs.reserve(_timer_N);
    // O(1) is easy to detect manually, plus it's difficult to detect it via
    // correlation, so just assume it's not that and move on
    ratio_type _timer_ns[_timer_N];
    // Logarithmic
    for (ratio_type _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
        _timer_ns[_timer_i] = std::ilogb(_timer_ns_orig[_timer_i]);
    }
    _timer_crs.emplace_back(std::make_pair(LOGARITHMIC,
            _timer_cc(_timer_times, _timer_ns)));
    // Linear
    for (ratio_type _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
        _timer_ns[_timer_i] = _timer_ns_orig[_timer_i];
    }
    _timer_crs.emplace_back(std::make_pair(LINEAR,
            _timer_cc(_timer_times, _timer_ns)));
    // Linearithmic
    for (ratio_type _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
        _timer_ns[_timer_i] = std::ilogb(_timer_ns_orig[_timer_i]) *
                _timer_ns_orig[_timer_i];
    }
    _timer_crs.emplace_back(std::make_pair(LINEARITHMIC,
            _timer_cc(_timer_times, _timer_ns)));
    // Quadratic
    for (ratio_type _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
        _timer_ns[_timer_i] = _timer_ns_orig[_timer_i] *
                _timer_ns_orig[_timer_i];
    }
    _timer_crs.emplace_back(std::make_pair(QUADRATIC,
            _timer_cc(_timer_times, _timer_ns)));
    // Cubic
    for (ratio_type _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
        ratio_type _timer_tmp = _timer_ns_orig[_timer_i];
        _timer_ns[_timer_i] = _timer_tmp * _timer_tmp * _timer_tmp;
    }
    _timer_crs.emplace_back(std::make_pair(CUBIC,
            _timer_cc(_timer_times, _timer_ns)));
    // Quartic
    for (ratio_type _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
        ratio_type _timer_tmp = _timer_ns_orig[_timer_i];
        _timer_ns[_timer_i] = (_timer_tmp * _timer_tmp) *
                (_timer_tmp * _timer_tmp);
    }
    _timer_crs.emplace_back(std::make_pair(QUADRATIC,
            _timer_cc(_timer_times, _timer_ns)));
    // Only calclate these next ones if they won't overflow
    // Exponential
    if (_timer_max_n < std::ilogb(std::numeric_limits<ratio_type>::max())) {
        for (ratio_type _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
            auto _timer_pow2 = [] (ratio_type _timer_e) {
                ratio_type _timer_result = 1;
                ratio_type _timer_b = 2;
                for ( ; _timer_e > 0; _timer_e /= 2) {
                    if (_timer_e % 2 == 1) {
                        _timer_result *= _timer_b;
                    }
                    _timer_b *= _timer_b;
                }
                return _timer_result;
            };
            _timer_ns[_timer_i] = _timer_pow2(_timer_ns_orig[_timer_i]);
        }
        _timer_crs.emplace_back(std::make_pair(EXPONENTIAL,
                _timer_cc(_timer_times, _timer_ns)));
    }
    // Factorial
    if (_timer_max_n <= 20) {
        for (ratio_type _timer_i = 0; _timer_i < _timer_N; ++_timer_i) {
            auto _timer_factorial = [] (ratio_type _timer_x) -> ratio_type {
                ratio_type _timer_f = 1;
                for (ratio_type _timer_j = 2; _timer_j <= _timer_x; ++_timer_j) {
                    _timer_f *= _timer_j;
                }
                return _timer_f;
            };
            _timer_ns[_timer_i] = _timer_factorial(_timer_ns_orig[_timer_i]);
        }
        _timer_crs.emplace_back(std::make_pair(FACTORIAL,
                _timer_cc(_timer_times, _timer_ns)));
    }
    auto _timer_iter = std::max_element(_timer_crs.begin(), _timer_crs.end(),
            [] (const std::pair<time_complexity, _timer_real> & _timer_elem1,
                const std::pair<time_complexity, _timer_real> & _timer_elem2)
            {
                return _timer_elem1.second < _timer_elem2.second;
            });
    return _timer_iter->first;
}

#endif // _timer_h_
