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
                particle.o cascade.o tree.o params.o main.o
_UKWO        = mainf.o upyth.o a.o siglookup.o addpart.o init.o output.o \
              string.o input.o tabinit.o paulibl.o angdis.o dectim.o \
     proppot.o ukw.o anndec.o delpart.o getmass.o upmerge.o detbal.o \
     getspin.o iso.o blockres.o ityp2pdg.o \
     dwidth.o jdecay2.o whichres.o boxprg.o hepchg.o \
     cascinit.o hepcmp.o make22.o hepnam.o  \
     coload.o numrec.o saveinfo.o scatter.o error.o
 
# VPATH = src:../UKW
HYDROO = $(patsubst %,$(ODIR)/%,$(_HYDROO))
UKWO = $(patsubst %,$(ODIR)/%,$(_UKWO))

TARGET = hadronSampler
#------------------------------------------------------------------------------

$(TARGET): $(UKWO) $(HYDROO)
		$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) -lgfortran
		@echo "$@ done"

clean:
		@rm -f $(ODIR)/*.o $(TARGET)

$(ODIR)/%.o: src.visc/%.cpp src.visc/const.h
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(ODIR)/%.o: ../UKW/%.f
		$(F) $(FFLAGS) -c $< -o $@
