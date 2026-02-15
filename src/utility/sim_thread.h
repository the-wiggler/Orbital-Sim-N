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
#include <unistd.h>
#ifdef __linux__
#include <sys/sysinfo.h>
#endif
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

// --- thread abstraction ---
typedef struct {
#ifdef _WIN32
    HANDLE handle;
#else
    pthread_t handle;
#endif
} thread_t;

#ifdef _WIN32
typedef DWORD (WINAPI *thread_func_t)(LPVOID);
#else
typedef void* (*thread_func_t)(void*);
#endif

static inline void thread_create(thread_t *thead, thread_func_t func, void *arg) {
#ifdef _WIN32
    thead->handle = CreateThread(NULL, 0, func, arg, 0, NULL);
#else
    pthread_create(&thead->handle, NULL, func, arg);
#endif
}

static inline void thread_join(thread_t *thread) {
#ifdef _WIN32
    WaitForSingleObject(thread->handle, INFINITE);
    CloseHandle(thread->handle);
#else
    pthread_join(thread->handle, NULL);
#endif
}

static inline int thread_get_cpu_count(void) {
#ifdef _WIN32
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return (int)si.dwNumberOfProcessors;
#elif defined(__linux__)
    return get_nprocs();
#elif defined(__APPLE__)
    int n = (int)sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? n : 4;
#else
    return 4;
#endif
}

#endif //ORBITSIMULATION_SIM_THEAD_H