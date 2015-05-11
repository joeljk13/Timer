#ifndef TIMER_H_
#define TIMER_H_ 1

/*
 * By default, the timer calls the function a few extra times that aren't
 * measured to get it into cache and ensure more consistent running times.
 * Specifying this option in flags will stop this.
 */
#define TIMER_NO_EXTRA 1

struct timer {
    unsigned long long ns;
    unsigned int reps;
};

/*
 * Measures function 'func'. Sets timer->ns to the number of nanoseconds it took,
 * and timer->reps to the number of repetitions. Uses the existing values of
 * timer->ns and timer->reps as maximums - it won't do any more repetitions or take
 * significantly more time than those specify. However, you can set one of them
 * to 0 to make it unlimited. 0 is returned on success, -1 on failure.
 */
int
timer_measure(void (*func)(void), struct timer *timer, int flags);

// These next functions are shortcuts that use timer_measure and use the
// default flags.

int
timer_measure_ms(void (*func)(void), unsigned long long ms,
                 struct timer *timer);

int
timer_measure_reps(void (*func)(void), unsigned int reps, struct timer *timer);

#endif
