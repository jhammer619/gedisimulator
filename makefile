# Makefile for iter_thresh
#
HOME =/Users/stevenhancock
LIBS = -lm -lgsl -lgslcblas -ltiff -lgeotiff
MINDIR=${HOME}/src/minpack
GSLDIR=${HOME}/src/GSL/gsl-1.16
INCLS = -I/usr/local/include -I${HOME}/src/headers -I$(MINDIR) -I${GSLDIR} -I${GSLDIR}/fft
CFLAGS += -Wall
#CFLAGS += -g
CFLAGS += -O3
LOCAL_OTHERS = libLasRead.o gaussFit.o libLasProcess.o tiffWrite.o
GSLFit=linear.o
MIN=mpfit.o

CC = gcc
#CC= /opt/SUNWspro/bin/cc

#CFLAGS += -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

THIS=gediRat

$(THIS):	$(THIS).o $(LOCAL_OTHERS) $(GSLDIR)/fit/$(GSLFIT) ${MINDIR}/$(MIN)
		$(CC) $(CFLAGS) $(LOCAL_OTHERS) $(GSLFIT) $(MIN) $@.o -o $@ $(LIBS) $(CFLAGS) $(INCLS)

.c.o:		$<
		$(CC) $(CFLAGS) -I. $(INCLS) -D$(ARCH)  -c $<

clean:
		\rm -f *% *~ *.o #$(TOOLDIR)/*.o

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/$(THIS).old
		cp $(THIS) $(HOME)/bin/$(ARCH)/$(THIS)


