#   ============================================================================
#   Youhua Xu Dec-13-2017
#   ============================================================================
#   Make a build dir for compilation
MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
BINDIR = $(MDIR)/bin

#	make build & binary dirs
.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; fi;
	touch $(WRKDIR)/.base;
	if ! [ -e $(BINDIR) ]; then mkdir $(BINDIR) ; fi;
	touch $(BINDIR)/.base;

#   ============================================================================
#   Set the source file path
vpath %.cpp main:source
vpath %.cc  main:source
vpath %.c   main:source
vpath %.hpp include
vpath %.h   include
vpath %.o build
vpath .base build

INCLUDES = -I $(MDIR)/include

CPP			= g++
CC          = mpic++
CCFLAG  	= -Wall -DHAVE_INLINE

OPTFLAG		= -O3 #-ffast-math #( not recognized by intel compiler )

LDFLAG      =
#   http://www.tuicool.com/articles/NBfeYj
ifeq ($(shell uname -s),Linux)
	LDFLAG	+= -Wl,--no-as-needed
endif
LIBS		= -limcmc
LIBS		+= -larmadillo -lblas -llapack
LIBS        += -lgsl -lgslcblas

ifeq ($(shell uname -s),Darwin)
	LDFLAG	+= -framework Accelerate #(-framework Accelerate is for Mac OSX)
endif

%.o: %.cpp .base
	cd $(WRKDIR);$(CC) $(OMPFLAG) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

%.o: %.cc .base
	cd $(WRKDIR);$(CC) $(OMPFLAG) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

%.o: %.c .base
	cd $(WRKDIR);$(CC) $(OMPFLAG) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

VERBOSE = Verbose.o

GSL     = GSL_Integration.o GSL_Interpolation.o GSL_GaussianRandom.o GSL_LinAlg.o

COSMO   = CosmoParams.o Cosmology.o BackgroundFunctions.o Background.o

SNE     = Cosmo_SNE.o Likelihood_SNE.o Data_SNE.o Data_SNE_Mock.o Data_SNE_UNION.o Data_SNE_LSST.o Data_SNE_SNLS.o Data_SNE_JLA.o ini.o
LIKE	= Likelihoods.o
ETA2E   = Likelihood_eta2E.o

MAIN_ISAMPLER	= iSampler.o

all:iSampler

#	iSampler replaces the old sampling exe 'iCosmo'
iSampler:${MAIN_ISAMPLER} ${COSMO} ${SNE} ${GSL} ${LIKE} ${ETA2E} ${VERBOSE}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

#   ================================================================================================
.PHONY:clean tidy run
clean: .base
	rm -rf $(WRKDIR);
tidy:
	make clean; rm -rf $(BINDIR);
