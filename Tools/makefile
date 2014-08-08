# Pull in ROOT
ROOT_INCLUDES=`root-config --cflags`
ROOT_ALL=`root-config --cflags --libs`
ROO_INCLUDES=${ROOT_ALL} -I${ROOFITSYS}/include -L${ROOFITSYS}/lib -lRooFit -lRooFitCore -lRooStats

#Compiler
CC=g++ -g -O2
CC11=${CC} -std=gnu++0x
CCC=${CC} -c

all: Joe

Joe: hnuPlots.C
	${CC} \
        ${ROO_INCLUDES} \
        -o hnuPlots.exe \
        hnuPlots.C

# All .o files rely on the .cpp files
%.o:%.cpp
	${CCC} ${ROO_INCLUDES} $< -o $@

# All .o files rely on the .h files
%.o:%.h
	${CCC} ${ROO_INCLUDES} $< -o $@

