#include "timer.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

static void
func1(void)
{
    printf("Hello world\n");
}

static void
func2(void)
{
    puts("Hello world");
}

int main(void)
{
    int err;
    unsigned int reps = 1024u * 1024u;
    unsigned long long ns = 30LL * 1000000000LL;
    struct timer timer1, timer2;

    timer1.ns = ns;
    timer1.reps = 0;
    timer2 = timer1;
    (void)reps;

    err = timer_measure(&func1, &timer1, 0);
    assert(err == 0);

    err = timer_measure(&func2, &timer2, 0);
    assert(err == 0);

    fprintf(stderr, "ns/rep = %llu, %llu\n", timer1.ns / timer1.reps, timer2.ns
        / timer2.reps);

    return 0;
}
