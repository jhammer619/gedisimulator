# Makefile for GEDI simulator tools

LIBS = -lm -lgsl -lgslcblas -ltiff -lgeotiff
INCLS = -I/usr/local/include -I$(HANCOCKTOOLS_ROOT) -I$(CMPFIT_ROOT) -I${GSL_ROOT} -I${GSL_ROOT}/fft -I${LIBCLIDAR_ROOT}
CFLAGS += -Wall
CFLAGS += -O3
LIBFILES = $(LIBCLIDAR_ROOT)/libLasProcess.o $(LIBCLIDAR_ROOT)/libLasRead.o $(LIBCLIDAR_ROOT)/tiffWrite.o $(LIBCLIDAR_ROOT)/gaussFit.o $(LIBCLIDAR_ROOT)/libLidVoxel.o  $(LIBCLIDAR_ROOT)/libTLSread.o
LOCLIB = libLasProcess.o libLasRead.o tiffWrite.o gaussFit.o libLidVoxel.o libTLSread.o
GSLFit=linear.o
MIN=mpfit.o
ARCH=$(shell uname -m)

CC = gcc

TARGET=gediRat

$(TARGET):	$(TARGET).o $(GSL_ROOT)/fit/$(GSLFIT) ${CMPFIT_ROOT}/$(MIN) ${LIBFILES}
		$(CC) $(CFLAGS) $(GSLFIT) ${CMPFIT_ROOT}/$(MIN) $(LOCLIB) $@.o -o $@ $(LIBS) $(CFLAGS) $(INCLS)

.c.o:		$<
		$(CC) $(CFLAGS) -I. $(INCLS) -D$(ARCH)  -c $<

clean:
		\rm -f *% *~ *.o

install:
		touch $(HOME)/bin/$(ARCH)/$(TARGET)
		mv $(HOME)/bin/$(ARCH)/$(TARGET) $(HOME)/bin/$(ARCH)/$(TARGET).old
		cp $(TARGET) $(HOME)/bin/$(ARCH)/$(TARGET)


