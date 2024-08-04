CPP = g++
CC = gcc
CFLAGS = -Wall $(shell pkg-config cfitsio --cflags) $(shell root-config --cflags) -O3
GCCNEWERTHAN47 := $(shell expr `gcc -dumpversion` \>= 4.7)
#add -std flag only if it doesn't already exist (on newer ROOT versions, root-config adds -std)
STD.0 = -std=c++0x
STD.1 = -std=c++11
ifeq (,$(findstring -std=,$(CFLAGS)))
    CFLAGS += $(STD.$(GCCNEWERTHAN47))
endif
LIBS = $(shell pkg-config cfitsio --libs) $(shell root-config --libs)
GLIBS = 
GLIBS +=
OBJECTS = skipper2root.o 
HEADERS = globalConstants.h

ALL : skipper2root.exe
	@echo "Listo!"

skipper2root.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o skipper2root.exe $(LIBS) $(GLIBS) $(CFLAGS)

skipper2root.o : skipper2root.cc $(HEADERS)
	$(CPP) -c skipper2root.cc -o skipper2root.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
