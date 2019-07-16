# compiler to use
CC = c++
CFLAGS=-O3 -std=c++14 
OFLAGS=-lpthread -fopenmp

LOCALPATH=/broad/IDP-Dx_work/nirmalya/local
INC = -I${LOCALPATH}/include

BOOSTLIBS=${LOCALPATH}/lib/libboost_serialization.a ${LOCALPATH}/lib/libboost_regex.a -lboost_iostreams -lz
PROG_OPT_LIB=${LOCALPATH}/lib/libboost_program_options.a 
STXXLLIBS=${LOCALPATH}/lib/libstxxl.a

all: tools
	
tools:
	$(CC) $(CFLAGS) $(INC) dual_bc_splitter.cpp -o dual_bc_splitter $(BOOSTLIBS) $(PROG_OPT_LIB) 

clean:
	rm -f *.o dual_bc_splitter
