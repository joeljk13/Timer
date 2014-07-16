#include <chrono>

#define BOOST
#include "timer.h"

#include <algorithm>
#include <iostream>
#include <deque>
#include <vector>
#include <random>

#include <cstdlib>
#include <ctime>

#define BUF_SIZE 50

// PR = Print Ratio
#define PR(x) \
do { \
    std::cout \
        << (#x) \
        << ": \n\tRatio: " \
        << (x).get_ratio() \
        << "\n\tCount: " \
        << (x).get_duration().count() \
        << std::endl; \
} while (0)

struct Type {
    Type() = default;
    Type(int i) : I (i) { }
    Type(const Type & t) = default;
    int I;
    char buf[BUF_SIZE];
    Type & operator =(const Type & t) = default;
    Type & operator =(int i) {
        I = i;
        return *this;
    }
    operator int() const {
        return I;
    }
    Type & operator +=(const Type & t) {
        I += t.I;
        return *this;
    }
    Type & operator +=(int i) {
        I += i;
        return *this;
    }
};

constexpr unsigned int REPS = 10;
constexpr std::chrono::seconds TM (15);

template <class T>
void test_push(std::deque<T> & deque, std::vector<T> & vector,
        std::vector<T> list)
{
    unsigned int num_elems = list.size();
    auto deque_push = [&] () {
        unsigned int n = num_elems;
        while (n--) {
            deque.push_back(list[n]);
        }
    };
    auto vector_push = [&] () {
        unsigned int n = num_elems;
        while (n--) {
            vector.push_back(list[n]);
        }
    };
    timer::timer deque_push_tm (deque_push, REPS);
    timer::timer vector_push_tm (vector_push, REPS);
    PR(deque_push_tm);
    PR(vector_push_tm);
}

template <class T>
void test_traverse(std::deque<T> & deque, std::vector<T> & vector)
{
    T n = 0;
    auto deque_traverse = [&] () {
        for (T i : deque) {
            n += i;
        }
    };
    auto vector_traverse = [&] () {
        for (T i : vector) {
            n += i;
        }
    };
    timer::timer deque_traverse_tm (deque_traverse, REPS);
    timer::timer vector_traverse_tm (vector_traverse, REPS);
    deque_traverse_tm.measure(deque_traverse, TM);
    vector_traverse_tm.measure(vector_traverse, TM);
    PR(deque_traverse_tm);
    PR(vector_traverse_tm);
    if (n == 0) {
        std::cout << std::endl;
    }
}

template <class T>
void test_at(std::deque<T> & deque, std::vector<T> & vector)
{
    T n = 0;
    auto deque_at = [&] () {
        unsigned int i = deque.size();
        while (i--) {
            n += deque.at(i);
        }
    };
    auto vector_at = [&] () {
        unsigned int i = vector.size();
        while (i--) {
            n += vector.at(i);
        }
    };
    timer::timer deque_at_tm (deque_at, REPS);
    timer::timer vector_at_tm (vector_at, REPS);
    deque_at_tm.measure(deque_at, TM);
    vector_at_tm.measure(vector_at, TM);
    PR(deque_at_tm);
    PR(vector_at_tm);
    if (n == 0) {
        std::cout << std::endl;
    }
}

template <class T>
void test_shuffle(std::deque<T> & deque, std::vector<T> & vector)
{
    auto deque_shuffle = [&] () {
        std::shuffle(deque.begin(), deque.end(),
                std::default_random_engine(std::rand()));
    };
    auto vector_shuffle = [&] () {
        std::shuffle(vector.begin(), vector.end(),
                std::default_random_engine(std::rand()));
    };
    timer::timer deque_shuffle_tm (deque_shuffle, REPS);
    timer::timer vector_shuffle_tm (vector_shuffle, REPS);
    deque_shuffle_tm.measure(deque_shuffle, TM);
    vector_shuffle_tm.measure(vector_shuffle, TM);
    PR(deque_shuffle_tm);
    PR(vector_shuffle_tm);
}

template <class T>
void test_sort(std::deque<T> & deque, std::vector<T> & vector)
{
    auto deque_sort = [&] () {
        std::deque<T> tmp (deque);
        std::sort(tmp.begin(), tmp.end());
    };
    auto vector_sort = [&] () {
        std::vector<T> tmp (vector);
        std::sort(tmp.begin(), tmp.end());
    };
    timer::timer deque_sort_tm (deque_sort, REPS);
    timer::timer vector_sort_tm (vector_sort, REPS);
    deque_sort_tm.measure(deque_sort, TM);
    vector_sort_tm.measure(vector_sort, TM);
    PR(deque_sort_tm);
    PR(vector_sort_tm);
}

template <class T>
void test(unsigned int num_elems)
{
    if (num_elems < 1000000) {
        return;
    }
    std::vector<T> list;
    list.reserve(num_elems);
    unsigned int num = num_elems;
    while (num--) {
        list.push_back(static_cast<T>(std::rand()));
    }
    std::deque<T> deque;
    std::vector<T> vector;
    std::cout << "Starting timing for " << num_elems << " elements" << std::endl;
    // This puts all the elements in, so it's necessary
    test_push(deque, vector, list);
    test_traverse(deque, vector);
    test_at(deque, vector);
    test_sort(deque, vector);
    test_shuffle(deque, vector);
}

int main()
{
    std::srand(std::time(NULL));
    std::cout << "USING unsigned char" << std::endl;
    for (int n = 1; n <= 1000000; n *= 10) {
        test<unsigned char>(n);
    }
    std::cout << "USING int" << std::endl;
    for (int n = 1; n <= 1000000; n*= 10) {
        test<int>(n);
    }
    std::cout << "USING unsigned int" << std::endl;
    for (int n = 1; n <= 1000000; n*= 10) {
        test<unsigned int>(n);
    }
    std::cout << "USING long long" << std::endl;
    for (int n = 1; n <= 1000000; n*= 10) {
        test<long long>(n);
    }
    std::cout << "USING double" << std::endl;
    for (int n = 1; n <= 1000000; n*= 10) {
        test<double>(n);
    }
    std::cout << "USING Type with buf = 50" << std::endl;
    for (int n = 1; n <= 1000000; n*= 10) {
        test<Type>(n);
    }
    return 0;
}

/*
void fn()
{
    std::cout << "This is a test message." << std::endl;
}

int main()
{
    const int bins = 50;
    const int count = 1000;
    const int repetitions = 1000;
    const char marker = '*';
    int **array = new int *[bins];
    for (int i = 0; i < bins; ++i) {
        array[i] = new int[count];
    }
    dist(&fn, count, array, bins, repetitions);
    std::cout << std::endl << "-----Actual Output-----" << std::endl <<
        std::endl;
    for (int i = 0; i < bins; ++i) {
        for (int j = 0; j < count; ++j) {
            if (array[i][j] > 0) {
                std::cout << marker;
            }
        }
        std::cout << std::endl << std::endl;
    }
    for (int i = 0; i < bins; ++i) {
        delete [] array[i];
    }
    delete [] array;
    return 0;
}
*/
