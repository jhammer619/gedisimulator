# Makefile for GEDI simulator tools

LIBS = -lm -lgsl -lgslcblas -ltiff -lgeotiff -L${GSL_ROOT}
INCLS = -I/usr/local/include -I$(HANCOCKTOOLS_ROOT) -I$(CMPFIT_ROOT) -I${LIBCLIDAR_ROOT}
CFLAGS += -Wall
CFLAGS += -O3
#CFLAGS += -g
LIBFILES = $(LIBCLIDAR_ROOT)/libLasProcess.o $(LIBCLIDAR_ROOT)/libLasRead.o $(LIBCLIDAR_ROOT)/tiffWrite.o $(LIBCLIDAR_ROOT)/gaussFit.o $(LIBCLIDAR_ROOT)/libLidVoxel.o  $(LIBCLIDAR_ROOT)/libTLSread.o
LOCLIB = libLasProcess.o libLasRead.o tiffWrite.o gaussFit.o libLidVoxel.o libTLSread.o
GSLFit=linear.o
MIN=mpfit.o
ARCH=$(shell uname -m)

CC = gcc

THIS=gediRat

$(THIS):	$(THIS).o ${CMPFIT_ROOT}/$(MIN) ${LIBFILES}
		$(CC) $(CFLAGS) $(MIN) $(LOCLIB) $@.o -o $@ $(LIBS) $(CFLAGS) $(INCLS)

.c.o:		$<
		$(CC) $(CFLAGS) -I. $(INCLS) -D$(ARCH)  -c $<

clean:
		\rm -f *% *~ *.o

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/$(THIS).old
		cp $(THIS) $(HOME)/bin/$(ARCH)/$(THIS)


