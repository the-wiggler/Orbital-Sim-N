//
// Created by java on 1/4/2026.
//

#ifndef ORBITSIMULATION_SIM_THEAD_H
#define ORBITSIMULATION_SIM_THEAD_H
#include "../globals.h"
#ifdef _WIN32
    #include <windows.h>
#else
#include <pthread.h>
#endif

typedef struct {
    union {
#ifdef _WIN32
        CRITICAL_SECTION win_cs;
#else
        pthread_mutex_t posix_mtx;
#endif
    } u;
} mutex_t;

// mutex lock function
static inline void mutex_lock(mutex_t *mutex) {
#ifdef _WIN32
    EnterCriticalSection(&mutex->u.win_cs);
#else
    pthread_mutex_lock(&mutex->u.posix_mtx);
#endif
}

// mutex unlock function
static inline void mutex_unlock(mutex_t *mutex) {
#ifdef _WIN32
    LeaveCriticalSection(&mutex->u.win_cs);
#else
    pthread_mutex_unlock(&mutex->u.posix_mtx);
#endif
}

// mutex initialization function
static inline void mutex_init(mutex_t *mutex) {
#ifdef _WIN32
    InitializeCriticalSection(&mutex->u.win_cs);
#else
    pthread_mutex_init(&mutex->u.posix_mtx, NULL);
#endif
}

// mutex destruction function
static inline void mutex_destroy(mutex_t *mutex) {
#ifdef _WIN32
    DeleteCriticalSection(&mutex->u.win_cs);
#else
    pthread_mutex_destroy(&mutex->u.posix_mtx);
#endif
}

#endif //ORBITSIMULATION_SIM_THEAD_H