#include "timer.h"

#include <assert.h>
#include <time.h>
#include <errno.h>

#include <stdio.h>

static int
measure_reps(void (*func)(void), struct timer *timer, int flags)
{
    unsigned int reps;
    clock_t begin, end;

    assert(timer->ns == 0);
    assert(timer->reps > 0);

    reps = timer->reps;

    if (!(flags & TIMER_NO_EXTRA)) {
        func();
    }

    begin = clock();
    do {
        func();
    }
    while (--reps);
    end = clock();

    assert(end >= begin);

    timer->ns = (unsigned long long)(end - begin) * (1000000000 /
        CLOCKS_PER_SEC);

    return 0;
}

static int
measure_ns_reps(void (*func)(void), struct timer *timer, int flags)
{
    // Find the minimum number of repetitions that measure to nonzero time
    // (call it n). Then use n to estimate how many repetitions will be
    // needed to pass half the given duration. After running the function
    // that many times, update the estimate of repitions for the next
    // quarter of 'duration'. Keep repeating this until the time alloted
    // has passed, or the number of repetitions estimated becomes less than
    // n.
    unsigned long long total_ns;
    clock_t start, stop, begin, end;
    unsigned int min_reps, total_reps;
    const unsigned int max_reps = timer->reps;
    const unsigned long long ns = timer->ns;

    assert(max_reps == 0);

    assert(timer->ns > 0);
    (void)flags;

    start = clock();
    stop = start + (clock_t)(ns / (1000000000 / CLOCKS_PER_SEC));
    begin = end = 0;
    total_reps = 0;
    total_ns = 0;
    min_reps = 1;

    for (;;) {
        unsigned int reps = min_reps;

        if (max_reps > 0 && total_reps + reps > max_reps) {
            reps = max_reps - total_reps;
        }

        begin = clock();
        do {
            func();
        }
        while (--reps);
        end = clock();

        total_reps += min_reps;

        if (begin == end) {
            continue;
        }

        total_ns += (unsigned long long)((end - begin) * (1000000000 / CLOCKS_PER_SEC));

        break;
    }

    if (max_reps > 0 && total_reps >= max_reps) {
        goto exit;
    }

    for (;;) {
        unsigned int reps;

        if (end >= stop) {
            break;
        }

        reps = (unsigned int)(((unsigned long long)((stop - end) * (1000000000
                        / CLOCKS_PER_SEC)) * total_reps) / (total_ns * 2));
        if (reps < min_reps) {
            break;
        }

        if (max_reps > 0 && total_reps + reps > max_reps) {
            if (total_reps >= max_reps) {
                break;
            }
            reps = max_reps - total_reps;
        }

        begin = clock();
        do {
            func();
        }
        while (--reps);
        end = clock();

        total_ns += (unsigned long long)((end - begin) / (1000000000 / CLOCKS_PER_SEC));
        total_reps += reps;
    }

exit:
    timer->ns = total_ns;
    timer->reps = total_reps;

    return 0;
}

int
timer_measure(void (*func)(void), struct timer *timer, int flags)
{
    if (func == NULL || timer == NULL) {
        errno = EINVAL;
        return -1;
    }

    if (timer->ns == 0) {
        return measure_reps(func, timer, flags);
    }

    return measure_ns_reps(func, timer, flags);
}

int
timer_measure_ms(void (*func)(void), unsigned long long ms,
                 struct timer *timer)
{
    timer->ns = ms * 1000;
    timer->reps = 0;

    if (timer_measure(func, timer, 0) != 0) {
        return -1;
    }

    return 0;
}

int
timer_measure_reps(void (*func)(void), unsigned int reps, struct timer *timer)
{
    timer->ns = 0;
    timer->reps = reps;

    if (timer_measure(func, timer, 0) != 0) {
        return -1;
    }

    return 0;
}
