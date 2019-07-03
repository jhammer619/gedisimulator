# Makefile for GEDI simulator tools

LIBS = -lm -lgsl -lgslcblas -ltiff -lgeotiff -lhdf5 -L${GSL_ROOT} -L${HDF5_LIB}/lib -lgdal
INCLS = -I/usr/local/include -I$(HANCOCKTOOLS_ROOT) -I$(CMPFIT_ROOT) -I${LIBCLIDAR_ROOT} -I${HDF5_LIB}/include -I. -I/usr/include/libgeotiff -I/usr/include/gdal
CFLAGS += -Wall 
#CFLAGS += -Wl,--verbose
#CFLAGS += -O3
CFLAGS += -g
LIBFILES = $(LIBCLIDAR_ROOT)/libLasProcess.o $(LIBCLIDAR_ROOT)/libLasRead.o $(LIBCLIDAR_ROOT)/tiffWrite.o $(LIBCLIDAR_ROOT)/gaussFit.o $(LIBCLIDAR_ROOT)/libLidVoxel.o  $(LIBCLIDAR_ROOT)/libTLSread.o  $(LIBCLIDAR_ROOT)/libLidarHDF.o $(GEDIRAT_ROOT)/gediIO.o $(LIBCLIDAR_ROOT)/libOctree.o $(GEDIRAT_ROOT)/gediNoise.o $(GEDIRAT_ROOT)/photonCount.o
LOCLIB = libLasProcess.o libLasRead.o tiffWrite.o gaussFit.o libLidVoxel.o libTLSread.o libLidarHDF.o gediIO.o libOctree.o gediNoise.o photonCount.o
GSLFit=linear.o
MIN=mpfit.o
ARCH=$(shell uname -m)

CC = gcc

THIS=gediRat

$(THIS):	$(THIS).o ${CMPFIT_ROOT}/$(MIN) ${LIBFILES}
		$(CC) $(CFLAGS) $(MIN) $(LOCLIB) $@.o -o $@ $(LIBS) $(INCLS)

.c.o:		$<
		$(CC) $(CFLAGS) -I. $(INCLS) -D$(ARCH)  -c $<

clean:
		\rm -f *% *~ *.o

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/$(THIS).old
		cp $(THIS) $(HOME)/bin/$(ARCH)/$(THIS)


