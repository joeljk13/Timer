#ifndef TIMER_H_
#define TIMER_H_ 1

/*
 * By default, the timer calls the function a few extra times that aren't
 * measured to get it into cache and ensure more consistent running times.
 * Specifying this option in flags will stop this.
 */
#define TIMER_NO_EXTRA 1

/*
 * By default, each of these functions assumes that the function is single
 * threaded. Specify this option to make sure that timing is done properly with
 * multi-threaded functions.
 */
#define TIMER_MULTI_THREAD 2

/*
 * The timer might create new processes to isolate the code being timed.
 * Specifying this flag prevents any new processes from being created.
 */
#define TIMER_NOPROC 4

/*
 * The timer might create new threads to isolate the code being timed.
 * Specifying this flag prevents any new threads from being created.
 */
#define TIMER_NOTHREAD 8

struct timer {
    unsigned long long ns;
    unsigned long reps;
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
// default flags. They just use 'timer' as an out argument.

int
timer_measure_ms(void (*func)(void), unsigned long long ms,
                 struct timer *timer);

int
timer_measure_reps(void (*func)(void), unsigned long reps, struct timer *timer);

#endif
