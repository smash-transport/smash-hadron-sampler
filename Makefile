ODIR           = obj

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = g++
F             = gfortran
CXXFLAGS      = -fPIC -O3
LD            = g++
LDFLAGS       = -O3
FFLAGS        = -fPIC $(ROOTCFLAGS) -O3

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

_HYDROO        = DecayChannel.o ParticlePDG2.o DatabasePDG2.o UKUtility.o gen.o \
                particle.o tree.o params.o main.o

HYDROO = $(patsubst %,$(ODIR)/%,$(_HYDROO))

TARGET = hadronSampler
#------------------------------------------------------------------------------

$(TARGET): $(HYDROO)
		$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) -lgfortran
		@echo "$@ done"

clean:
		@rm -f $(ODIR)/*.o $(TARGET)

$(ODIR)/%.o: src/%.cpp src/const.h
		$(CXX) $(CXXFLAGS) -c $< -o $@
