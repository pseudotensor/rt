CC=g++
CXX=g++
CFLAGS=-fopenmp -Ofast
CXXFLAGS=-fopenmp -Ofast # -pg -fprofile-arcs -ftest-coverage
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
