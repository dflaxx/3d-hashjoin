// from Guido (/home/src/moer/infra/)

#ifdef __cplusplus
extern "C" {
#endif

#ifndef INFRA_CBIND_TO_HW_THREAD
#define INFRA_CBIND_TO_HW_THREAD

// only enabled if not on Mac
#ifndef __MACH__
  // _GNU_SOURCE: https://stackoverflow.com/questions/5582211/what-does-define-gnu-source-imply
  #ifndef _GNU_SOURCE
    #define _GNU_SOURCE
  #endif
  #include <sched.h>
#endif

int cbind_to_hw_thread(const int aHwThreadNo, const int aMsg);
int cbind_sw_to_hw_thread(const pid_t aSwThreadNo, const int aHwThreadNo, const int aMsg);

#endif  // #ifndef INFRA_CBIND_TO_HW_THREAD

#ifdef __cplusplus
}
#endif
