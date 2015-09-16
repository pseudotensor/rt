CC=g++
CXX=g++
WARP=fast # can be: 1,2,3,fast # SEG-FAULTS with 0 :-S
# CFLAGS=-fopenmp -Ofast
# CXXFLAGS=-fopenmp -Ofast # -pg -fprofile-arcs -ftest-coverage
CFLAGS=-fopenmp -O${WARP}
CXXFLAGS=-fopenmp -O${WARP} # -pg -fprofile-arcs -ftest-coverage
# CFLAGS=-fopenmp 
# CXXFLAGS=-fopenmp # -pg -fprofile-arcs -ftest-coverage
LIBS=
RM=/bin/rm

EXECS = ASTRORAY_main
OBJS = ASTRORAY_main.o

$(EXECS): $(OBJS) Makefile

	$(CC) $(CFLAGS) -o $(EXECS) $(OBJS) $(LIBS)
$(OBJS): Makefile
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c
clean:
	$(RM) -f $(EXECS) $(OBJS) gmon.out *~ *.o *.oo
