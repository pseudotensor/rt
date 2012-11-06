CC=icc
CXX=icc
CFLAGS=-openmp
CXXFLAGS=-openmp
LIBS=-lstdc++
RM=/bin/rm

EXECS = transfer_quad
OBJS = transfer_quad.o

$(EXECS): $(OBJS) Makefile

	$(CC) $(CFLAGS) -o $(EXECS) $(OBJS) $(LIBS)
$(OBJS): Makefile
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c
clean:
	$(RM) -f $(EXECS) $(OBJS) gmon.out *~ *.o *.oo
                        