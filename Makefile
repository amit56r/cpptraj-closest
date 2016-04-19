# Makefile for closest
all: config.h install

install: a.out

config.h:
	@echo "config.h not present. Please run ./configure"
	@exit 1

include config.h

SOURCES=Atom.cpp Box.cpp CpptrajStdio.cpp MaskToken.cpp NameType.cpp Vec3.cpp \
        AtomMask.cpp CoordinateInfo.cpp Frame.cpp Matrix_3x3.cpp  Residue.cpp \
        main.cpp ArgList.cpp StringRoutines.cpp DistRoutines.cpp Topology.cpp \
        CharMask.cpp FileName.cpp Range.cpp Action_Closest.cpp \
        Parm_Amber.cpp TextFormat.cpp CpptrajFile.cpp FileIO_Bzip2.cpp \
        FileIO_Gzip.cpp FileIO_Std.cpp NetcdfFile.cpp Traj_AmberNetcdf.cpp \
        cuda_kernels/core_kernels.cu \
        cuda_kernels/kernel_wrappers.cu \
        cuda_kernels/setup_routines.cpp 








OBJECTS=$(SOURCES:.cpp=.o)

a.out: $(OBJECTS)
	$(CXX) -o a.out $(OBJECTS) $(NETCDF_LIB_DIR) $(LDFLAGS)

.cpp.o:
	$(CXX) $(DEFINES) $(CXXFLAGS) -c -o $@ $<

NetcdfFile.o: NetcdfFile.cpp
	$(CXX) $(DEFINES) $(NETCDF_INC_DIR) $(CXXFLAGS) -c -o NetcdfFile.o NetcdfFile.cpp

Traj_AmberNetcdf.o: Traj_AmberNetcdf.cpp
	$(CXX) $(DEFINES) $(NETCDF_INC_DIR) $(CXXFLAGS) -c -o Traj_AmberNetcdf.o Traj_AmberNetcdf.cpp 

findDepend: FindDepend.o
	$(CXX) -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCES) > cpptrajdepend

dependclean:
	/bin/rm -f FindDepend.o findDepend

clean:
	/bin/rm -f $(OBJECTS) a.out

uninstall: clean dependclean
	/bin/rm -f config.h
	/bin/rm -f test.out

test: a.out
	./a.out > test.out
	diff test.out.save test.out

include cpptrajdepend
