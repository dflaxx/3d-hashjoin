HOSTNAME = $(shell cat /etc/hostname)

DO_OPTIMIZE = 1
DO_DEBUGINFO = 0
DO_ASSERTS = 0
DO_PROFILE_GPROF = 0
DO_SANITIZERS = 0

$(info ******* BUILD CONFIG ********)
$(info OPT  DEBUG  ASSERT  PROF  SAN)
$(info ${DO_OPTIMIZE}    ${DO_DEBUGINFO}      ${DO_ASSERTS}       ${DO_PROFILE_GPROF}     ${DO_SANITIZERS}) 
$(info host: ${HOSTNAME})
$(info *****************************)

CC       := clang
CXX      := clang++
#CC       := gcc
#CXX      := g++

CXXSTD = -std=c++20

ifeq (1, ${DO_OPTIMIZE})
CFLAGS := -Wall -Wextra -O3
CXXFLAGS := $(CXXSTD) -Wall -Wextra -Weffc++ -O3
else
$(warning WARNING: Optimization disabled!)
CFLAGS := -Wall -Wextra -O0
CXXFLAGS := $(CXXSTD) -Wall -Wextra -Weffc++ -O0
endif

ifeq (1, ${DO_DEBUGINFO})
$(warning WARNING: Compiling with debug info: -g3!)
CFLAGS := ${CFLAGS} -g3
CXXFLAGS := ${CXXFLAGS} -g3
else
CFLAGS := ${CFLAGS} -g0
CXXFLAGS := ${CXXFLAGS} -g0
endif

ifeq (1, ${DO_ASSERTS})
CFLAGS := ${CFLAGS} -DDEBUG
CXXFLAGS := ${CXXFLAGS} -DDEBUG
else
$(warning WARNING: Asserts are disabled!)
CFLAGS := ${CFLAGS} -DNDEBUG
CXXFLAGS := ${CXXFLAGS} -DNDEBUG
endif

ifeq (1, ${DO_PROFILE_GPROF})
CFLAGS := ${CFLAGS} -pg
CXXFLAGS := ${CXXFLAGS} -pg
endif

DEPFLAGS = -MMD -MP

ifeq (1, ${DO_SANITIZERS})
SANITIZERS=-fsanitize=address,undefined
endif

CTAGS      = ctags
TAGFLAGS   = --c-kinds=+p --fields=+iaS --extras=+q
# tag all .cc and .hh files, but do not expand directory .git
TAGFILES   = $(shell find . -type f -regex '.*\.\(hh\|cc\|h\|c\)' -not -path './.git/*')

TARGETS  = main_experiment1.out \
	   main_experiment4.out \
	   main_algebra_example.out
MAINOBJS = $(TARGETS:.out=.o)
UTILOBJS = util/GenRandIntVec.o util/csv_writer.o util/cbind_to_hw_thread.o
OBJS     = $(UTILOBJS) ht_statistics.o 

.PHONY: all clean printenv

all: $(TARGETS) tags

# Linking
$(TARGETS): %.out: %.o $(OBJS)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) $(SANITIZERS) $^ -o $@ 

# Compilation
%.o: %.cc
	@echo "Compiling $@ (C++ compiler)"
	$(CXX) $(CXXFLAGS) $(SANITIZERS) $(DEPFLAGS) $(IPATHS) -c -o $@ $<

%.o: %.c
	@echo "Compiling $@ (C compiler)"
	$(CC) $(CFLAGS) $(DEPFLAGS) $(IPATHS) -c -o $@ $<

-include *.d

clean:
	-rm -f $(TARGETS) $(MAINOBJS) $(OBJS) tags
	-find . -type f -name "*.d" -delete

tags: $(TAGFILES)
	@echo "Generating ctags..."
	-$(CTAGS) $(TAGFLAGS) $(TAGFILES)

printenv:
	@echo "HOSTNAME         = ${HOSTNAME}"
	@echo "DO_OPTIMIZE      = ${DO_OPTIMIZE}"
	@echo "DO_DEBUGINFO     = ${DO_DEBUGINFO}"
	@echo "DO_ASSERTS       = ${DO_ASSERTS}"
	@echo "DO_PROFILE_GPROF = ${DO_PROFILE_GPROF}"
	@echo "DO_SANITIZERS    = ${DO_SANITIZERS}"
	@echo "CC               = $(CC)"
	@echo "CFLAGS           = $(CFLAGS)"
	@echo "CXX              = $(CXX)"
	@echo "CXXFLAGS         = $(CXXFLAGS)"
	@echo "DEPFLAGS         = $(DEPFLAGS)"
	@echo "TARGETS          = $(TARGETS)"
	@echo "MAINOBJS         = $(MAINOBJS)"
	@echo "OBJS             = $(OBJS)"
	@echo "CTAGS            = $(CTAGS)"
	@echo "TAGFLAGS         = $(TAGFLAGS)"
	@echo "TAGFILES         = $(TAGFILES)"

