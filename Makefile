SHELL=/bin/sh

SRCS= \
GF.cpp \
Halo.cpp \
LSS.cpp \
NFW.cpp \
PS.cpp \
TF.cpp \
MF.cpp \
PDF.cpp \
CDF.cpp \
grid.cpp

HDRS= \
GF.h \
Halo.h \
LSS.h \
NFW.h \
PS.h \
TF.h \
MF.h \
PDF.h \
CDF.h \
grid.h

UTILS= pdf nongaussian convolve integrate integrate_std integrate_weighted difference test_halo test_mf test_ps test_ps_var test_all

###########################################################################
# Commands and options for compiling
########################################################################### 
OBJS= $(addsuffix .o, $(basename $(SRCS)))
 
CC= g++
CFLAGS=  -g
#WARNFLAGS= -Wall -W -Wshadow -fno-common
#MOREFLAGS= -ansi -pedantic -Wpointer-arith -Wcast-qual -Wcast-align \
#           -Wwrite-strings -fshort-enums 
LDFLAGS= -lgsl -lgslcblas -lgomp

###########################################################################
# Instructions to compile and link -- allow for different dependencies
########################################################################### 
 
%.o: %.cpp $(HDRS) $(MAKEFILE)
	$(CC) $(CFLAGS) $(WARNFLAGS) -c $< -o $@
    
pdf: main.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o pdf $(OBJS) main.cpp $(LDFLAGS) $(LIBS)

nongaussian: nongaussian_main.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o nongaussian $(OBJS) nongaussian_main.cpp $(LDFLAGS) $(LIBS)
	
convolve: convolve.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o convolve convolve.cpp $(CFLAGS) $(LDFLAGS) $(LIBS) 
	
integrate: integrate.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o integrate $(OBJS) integrate.cpp $(CFLAGS) $(LDFLAGS) $(LIBS) 

integrate_weighted: integrate_weighted.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o integrate_weighted $(OBJS) integrate_weighted.cpp $(CFLAGS) $(LDFLAGS) $(LIBS) 

integrate_std: integrate_std.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o integrate_std $(OBJS) integrate_std.cpp $(CFLAGS) $(LDFLAGS) $(LIBS) 
	
difference: difference.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o difference $(OBJS) difference.cpp $(CFLAGS) $(LDFLAGS) $(LIBS) 

test_halo: test_halo.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o test_halo $(OBJS) test_halo.cpp $(CFLAGS) $(LDFLAGS) $(LIBS)

test_mf: test_mf.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o test_mf $(OBJS) test_mf.cpp $(CFLAGS) $(LDFLAGS) $(LIBS)

test_ps: test_ps.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o test_ps $(OBJS) test_ps.cpp $(CFLAGS) $(LDFLAGS) $(LIBS)

test_ps_var: test_ps_var.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o test_ps_var $(OBJS) test_ps_var.cpp $(CFLAGS) $(LDFLAGS) $(LIBS)

test_all: test_all.cpp $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o test_all $(OBJS) test_all.cpp $(CFLAGS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(UTILS)
