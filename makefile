# Makefile for GEDI simulator tools

LIBS = -L${GSL_ROOT}/lib/$(ARCH) -L${HDF5_LIB}/lib -L${GDAL_LIB} #-L/anaconda3/lib

ifeq ($(OS),Windows_NT)
    CFLAGS += -DWIN32 -D_WIN32 -D_CRT_SECURE_NO_WARNINGS
	LIBS += -lgdal_i -lgsl -lgslcblas -ltiff -lgeotiff_i -lhdf5 -lhdf5_hl
else
	LIBS += -lgdal -lm -lgsl -lgslcblas -ltiff -lgeotiff -lhdf5
endif

INCLS = -I/usr/local/include -I$(HANCOCKTOOLS_ROOT) -I$(CMPFIT_ROOT) -I${LIBCLIDAR_ROOT} -I. -I${LIB_GEOTIFF} -I/usr/include/gdal -I${GSL_ROOT}/include  -I${HDF5_LIB}/include #-I/anaconda3/include
CFLAGS += -Wall -DUSEPHOTON
#CFLAGS += -Wl,--verbose
CFLAGS += -O3 -DM_PI=3.14159265358979323846
#CFLAGS += -g
LIBFILES = $(HANCOCKTOOLS_ROOT)/tools.o $(HANCOCKTOOLS_ROOT)/functionWrappers.o $(LIBCLIDAR_ROOT)/libLasProcess.o $(LIBCLIDAR_ROOT)/libLasRead.o $(LIBCLIDAR_ROOT)/tiffWrite.o $(LIBCLIDAR_ROOT)/gaussFit.o $(LIBCLIDAR_ROOT)/libLidVoxel.o  $(LIBCLIDAR_ROOT)/libTLSread.o  $(LIBCLIDAR_ROOT)/libLidarHDF.o gediIO.o $(LIBCLIDAR_ROOT)/libOctree.o gediNoise.o photonCount.o
LOCLIB = libLasProcess.o libLasRead.o tiffWrite.o gaussFit.o libLidVoxel.o libTLSread.o libLidarHDF.o gediIO.o libOctree.o gediNoise.o photonCount.o
GSLFit=linear.o
MIN=mpfit.o


CC = gcc

ifndef THIS
	THIS=gediRat
endif

$(THIS): ${CMPFIT_ROOT}/$(MIN) ${LIBFILES} $(THIS).o
		$(CC) $(CFLAGS) $(INCLS) $^ -o $@ $(LIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLS) $< -o $@

clean:
		\rm -f *% *~ $(LIBFILES)

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/
