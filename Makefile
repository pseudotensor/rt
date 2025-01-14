CC=icc
CXX=icc
CFLAGS=-openmp
CXXFLAGS=-openmp
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
