// from Guido (/home/src/moer/infra/)

// _GNU_SOURCE: https://stackoverflow.com/questions/5582211/what-does-define-gnu-source-imply
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "cbind_to_hw_thread.h"

#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <assert.h>

// only enabled on Mac
#ifdef __MACH__

/*
 * on MacOS
 * sysctlbyname(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen);
 * can be used to derive some system information
 * command:
 * sysctl -a 
 */

#include <sys/types.h>
#include <sys/sysctl.h>

#define SYSCTL_CORE_COUNT   "machdep.cpu.core_count"

typedef struct cpu_set {
  uint32_t    count;
} cpu_set_t;

static inline void
CPU_ZERO(cpu_set_t *cs) {
  cs->count = 0;
}

static inline void
CPU_SET(int num, cpu_set_t *cs) {
  cs->count |= (1 << num);
}

static inline int
CPU_ISSET(int num, cpu_set_t *cs) {
  return (cs->count & (1 << num));
}

int sched_getaffinity(pid_t pid, size_t cpu_size, cpu_set_t *cpu_set);
int sched_setaffinity(int aPid, size_t aSize, cpu_set_t* aMask);


int
sched_getaffinity(pid_t pid, size_t cpu_size, cpu_set_t *cpu_set) {
  int32_t core_count = 0;
  size_t  len = sizeof(core_count);
  int ret = sysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
  if (ret) {
    printf("error while get core count %d\n", ret);
    return -1;
  }
  cpu_set->count = 0;
  for (int i = 0; i < core_count; i++) {
    cpu_set->count |= (1 << i);
  }

  return 0;
}

int 
sched_setaffinity(int aPid, size_t aSize, cpu_set_t* aMask) {
  return 0;
}

/*
  Q: I don't think that's possible to link a thread with a specific core with Mac OS X …
  A: That’s correct.
  Q: To test, I've wrote a simple program to validate and 
     it seems impossible to print from which core the thread is executing.
  A: That’s also correct.
  A: Well, there are ways to work this out but the previous point means 
     it doesn’t do you any good (by the time you look at the result, you might be running on a different CPU).

   THREAD_AFFINITY_POLICY
   is pretty well covered in the header comments within
   <Kernel/thread_policy.h>
   To wit: This policy is experimental.

   citation from <Kernel/thread_policy.h>
   This may be used to express affinity relationships between threads in
   the task. Threads with the same affinity tag will be scheduled to
   share an L2 cache if possible. That is, affinity tags are a hint to
   the scheduler for thread placement.
   
   The namespace of affinity tags is generally local to one task.
   However, a child task created after the assignment of affinity tags by
   its parent will share that namespace. In particular, a family of
   forked processes may be created with a shared affinity namespace.
*/

#endif  // #ifdef __MACH__


int
cbind_to_hw_thread(const int aHwThreadNo, const int aMsg) {
  cpu_set_t lMask;
  CPU_ZERO(&lMask);
  CPU_SET(aHwThreadNo, &lMask);
  int lRc = sched_setaffinity(0, sizeof(lMask), &lMask);
  if(0 == lRc) { 
    if(0 < aMsg) {
      printf("# binding thread to hw-thread %d succeeded.\n", aHwThreadNo);
    }
  } else {
    if(aMsg) {
      printf("# binding thread to hw-thread %d failed: return code: %d.\n", aHwThreadNo, lRc);
    } else
    if(0 > aMsg) {
      assert(0 == aMsg);
    }
  }
  return (0 == lRc);
}

int
cbind_sw_to_hw_thread(const pid_t aSwThreadNo, const int aHwThreadNo, const int aMsg) {
  cpu_set_t lMask;
  CPU_ZERO(&lMask);
  CPU_SET(aHwThreadNo, &lMask);
  int lRc = sched_setaffinity(aSwThreadNo, sizeof(lMask), &lMask);
  if(0 == lRc) {
    if(0 < aMsg) {
      printf("# binding thread to hw-thread %d succeeded.\n", aHwThreadNo);
    }
  } else {
    if(aMsg) {
      printf("# binding thread to hw-thread %d failed: return code: %d.\n", aHwThreadNo, lRc);
    } else
    if(0 > aMsg) {
      assert(0 == aMsg);
    }
  }
  return (0 == lRc);
}



