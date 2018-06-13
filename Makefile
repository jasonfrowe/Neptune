#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F90 = ifort
F77 = ifort
#compiling object file flags
OPT1 = -O3 -mkl
#OPT1 = -O0 -g -CB -openmp -mkl
#OPT1 = -O3
OPT2 = -heap-arrays 0 
#OPT2 = 
FFLAGS = -c $(OPT1) $(OPT2)
#FFLAGS = -c -O0 -g -CB -warn $(OPT2)
#linking flags
LFLAGS = $(OPT1) $(OPT2)
#LFLAGS = -O0 -g -CB -warn $(OPT2)
#testing flags
#LFLAGS = -O0 -g -CB
#fitsio libraries
FITSIODIR = /usr/local/lib
#Pgplot plot libraries
PGPLOTDIR = /usr/local/lib
#X11 libraries
X11DIR = /usr/X11R6/lib
# libraries for linking PGPLOT
PLPLOTDIR = -I/usr/local/Cellar/plplot/5.9.11/lib/fortran/modules/plplot -I/usr/local/Cellar/plplot/5.9.11/include/plplot -L/usr/local/Cellar/plplot/5.9.11/lib
#LIBS = $(PLPLOTDIR) -lplplotf95d -lplplotf95cd
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/bin/
#utils source directory
UTILS = utils/

#Listing of programs to create.
all: fitdatav2 pixelfit fftpow

fftplotincl = 
fftplot: fftplot.f90 $(fftplotinc)
	$(F90) $(LFLAGS) -o $(BIN)$@ fftplot.f90 $(fftplotincl) $(LIBS2)

figure3incl = 
figure3: figure3.f90 $(figure3incl)
	$(F90) $(LFLAGS) -o $(BIN)$@ figure3.f90 $(figure3incl) $(LIBS2)

joinpartsincl = precision.o fitneptunepos.o lfit.o
joinparts: joinparts.f90 $(joinpartsincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ joinparts.f90 $(joinpartsincl)

gendataincl = precision.o ran2.o
gendata: gendata.f90 $(gendataincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ gendata.f90 $(gendataincl)

fftpowincl = precision.o fittingfft.o readfftdata.o plotdatascatter.o fftspec.o ran2.o rqsort.o makekernel.o resample.o plotspec.o fftstats.o stdev.o poorwavelet.o heatlut.o ovrwrt.o fitfft.o minpack.o 
fftpow: fftpow.f90 $(fftpowincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ fftpow.f90 $(fftpowincl) $(LIBS) -lfftw3

pixelfitincl = precision.o getdata.o plotdatascatter.o makekernel.o plotsamples.o fitneptunepos.o lfit.o pshapemodel.o pfitter.o lbfgsb.o timer.o linpack.o
pixelfit: pixelfit.f90 $(pixelfitincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ pixelfit.f90 $(pixelfitincl) $(LIBS)

fitdatav3incl = precision.o getdata.o cutoutliers.o stdev.o meddiff.o rqsort.o
fitdatav3: fitdatav3.f90 $(fitdatav3incl)
	$(F90) $(LFLAGS) -o $(BIN)$@ fitdatav3.f90 $(fitdatav3incl)

fitdatav2incl = precision.o getdata.o findjumps.o makekernel.o cutoutliers.o stdev.o fitline.o fitterv2.o pixelmodelv2.o lbfgsb.o timer.o gradient.o spcor.o spline.o exportdata.o linpack.o meddiff.o rqsort.o
fitdatav2: fitdatav2.f90 $(fitdatav2incl)
	$(F90) $(LFLAGS) -o $(BIN)$@ fitdatav2.f90 $(fitdatav2incl)

fitdataincl = precision.o fittingmod.o getdata.o plotdata.o fitline.o plotline.o makekernel.o displaykernel.o heatlut.o stdev.o rqsort.o lapack.o blas.o plotsamples.o plotdatascatter.o fitter.o lfit.o fitneptunepos.o minpack.o fcn.o pixelmodel.o cutoutliers.o
fitdata: fitdata.f90 $(fitdataincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ fitdata.f90 $(fitdataincl) $(LIBS) 

#building object libraries
meddiff.o: $(UTILS)meddiff.f90
	$(F90) $(FFLAGS) $(UTILS)meddiff.f90
fittingfft.o: $(UTILS)fittingfft.f90
	$(F90) $(FFLAGS) $(UTILS)fittingfft.f90
fitfft.o: $(UTILS)fitfft.f90
	$(F90) $(FFLAGS) $(UTILS)fitfft.f90
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
getdata.o: $(UTILS)getdata.f90
	$(F90) $(FFLAGS) $(UTILS)getdata.f90
plotdata.o: $(UTILS)plotdata.f90
	$(F90) $(FFLAGS) $(UTILS)plotdata.f90
plotdatascatter.o: $(UTILS)plotdatascatter.f90
	$(F90) $(FFLAGS) $(UTILS)plotdatascatter.f90
fitline.o: $(UTILS)fitline.f90
	$(F90) $(FFLAGS) $(UTILS)fitline.f90
plotline.o: $(UTILS)plotline.f90
	$(F90) $(FFLAGS) $(UTILS)plotline.f90
makekernel.o: $(UTILS)makekernel.f90
	$(F90) $(FFLAGS) $(UTILS)makekernel.f90
displaykernel.o: $(UTILS)displaykernel.f90
	$(F90) $(FFLAGS) $(UTILS)displaykernel.f90
heatlut.o: $(UTILS)heatlut.f90
	$(F90) $(FFLAGS) $(UTILS)heatlut.f90
stdev.o: $(UTILS)stdev.f
	$(F90) $(FFLAGS) $(UTILS)stdev.f
rqsort.o: $(UTILS)rqsort.f
	$(F90) $(FFLAGS) $(UTILS)rqsort.f
lapack.o: $(UTILS)lapack.f
	$(F90) $(FFLAGS) $(UTILS)lapack.f
blas.o: $(UTILS)blas.f
	$(F90) $(FFLAGS) $(UTILS)blas.f
plotsamples.o: $(UTILS)plotsamples.f90
	$(F90) $(FFLAGS) $(UTILS)plotsamples.f90
fitter.o: $(UTILS)fitter.f90
	$(F90) $(FFLAGS) $(UTILS)fitter.f90
lfit.o: $(UTILS)lfit.f
	$(F90) $(FFLAGS) $(UTILS)lfit.f
fitneptunepos.o: $(UTILS)fitneptunepos.f90
	$(F90) $(FFLAGS) $(UTILS)fitneptunepos.f90
minpack.o: $(UTILS)minpack.f
	$(F90) $(FFLAGS) $(UTILS)minpack.f
fcn.o: $(UTILS)fcn.f90
	$(F90) $(FFLAGS) $(UTILS)fcn.f90
fittingmod.o: $(UTILS)fittingmod.f90
	$(F90) $(FFLAGS) $(UTILS)fittingmod.f90
pixelmodel.o: $(UTILS)pixelmodel.f90
	$(F90) $(FFLAGS) $(UTILS)pixelmodel.f90
findjumps.o: $(UTILS)findjumps.f90
	$(F90) $(FFLAGS) $(UTILS)findjumps.f90
cutoutliers.o: $(UTILS)cutoutliers.f90
	$(F90) $(FFLAGS) $(UTILS)cutoutliers.f90
medfit.o: $(UTILS)medfit.f
	$(F90) $(FFLAGS) $(UTILS)medfit.f
fitterv2.o: $(UTILS)fitterv2.f90
	$(F90) $(FFLAGS) $(UTILS)fitterv2.f90
amoeba.o: $(UTILS)amoeba.f
	$(F90) $(FFLAGS) $(UTILS)amoeba.f
pixelmodelv2.o: $(UTILS)pixelmodelv2.f90
	$(F90) $(FFLAGS) $(UTILS)pixelmodelv2.f90
fittingmodv2.o: $(UTILS)fittingmodv2.f90
	$(F90) $(FFLAGS) $(UTILS)fittingmodv2.f90
lbfgsb.o: $(UTILS)lbfgsb.f
	$(F90) $(FFLAGS) $(UTILS)lbfgsb.f
linpack.o: $(UTILS)linpack.f
	$(F90) $(FFLAGS) $(UTILS)linpack.f
timer.o: $(UTILS)timer.f
	$(F90) $(FFLAGS) $(UTILS)timer.f
gradient.o: $(UTILS)gradient.f90
	$(F90) $(FFLAGS) $(UTILS)gradient.f90
spline.o: $(UTILS)spline.f
	$(F90) $(FFLAGS) $(UTILS)spline.f
spcor.o: $(UTILS)spcor.f90
	$(F90) $(FFLAGS) $(UTILS)spcor.f90
exportdata.o: $(UTILS)exportdata.f90
	$(F90) $(FFLAGS) $(UTILS)exportdata.f90
readfftdata.o: $(UTILS)readfftdata.f90
	$(F90) $(FFLAGS) $(UTILS)readfftdata.f90
fftspec.o: $(UTILS)fftspec.f90
	$(F90) $(FFLAGS) $(UTILS)fftspec.f90
ran2.o: $(UTILS)ran2.f
	$(F90) $(FFLAGS) $(UTILS)ran2.f
resample.o: $(UTILS)resample.f90
	$(F90) $(FFLAGS) $(UTILS)resample.f90
plotspec.o: $(UTILS)plotspec.f90
	$(F90) $(FFLAGS) $(UTILS)plotspec.f90
fftstats.o: $(UTILS)fftstats.f90
	$(F90) $(FFLAGS) $(UTILS)fftstats.f90
poorwavelet.o: $(UTILS)poorwavelet.f90
	$(F90) $(FFLAGS) $(UTILS)poorwavelet.f90
ovrwrt.o: $(UTILS)ovrwrt.f
	$(F90) $(FFLAGS) $(UTILS)ovrwrt.f
pshapemodel.o: $(UTILS)pshapemodel.f90
	$(F90) $(FFLAGS) $(UTILS)pshapemodel.f90
pfitter.o: $(UTILS)pfitter.f90
	$(F90) $(FFLAGS) $(UTILS)pfitter.f90

# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
