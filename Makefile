
PROGRAM = main.exe

ifndef WINE
CC      = gcc
WRES    = windres
else
CC      = winegcc
WRES    = wrc
endif

ifdef DEBUG
ODIR    = debug
CFLAGS  = -g -D_DEBUG -mno-cygwin
else
ODIR    = release
CFLAGS  = -Wall
endif

SRCDIR  = .
LIBDIRS =
INCDIRS =
LIBS    = -lm

SRCS = $(SRCDIR)/main.c\
	   $(SRCDIR)/cspline.c


OBJS = $(ODIR)/main.obj\
	   $(ODIR)/cspline.obj


HDRS = $(SRCDIR)/cspline.h


ALLOBJS  = $(OBJS)
ALLBIN   = $(ODIR)/$(PROGRAM)

all: $(ODIR)/$(PROGRAM)

cleanobjs:
	rm -f $(ALLOBJS)

cleanbin:
	rm -f $(ALLBIN)

clean: cleanobjs cleanbin

cleanall: cleanobjs cleanbin

$(ODIR)/$(PROGRAM): $(OBJS) $(HDRS)
	$(CC) -o $(ODIR)/$(PROGRAM) $(OBJS) $(INCDIRS) $(LIBDIRS) $(LIBS)

$(ODIR)/%.res : %.rc $(HDRS) $(ODIR)
	$(WRES) --use-temp-file -O coff $< $@

$(ODIR)/%.obj : %.cpp $(HDRS) $(ODIR)
	$(CC) $(CFLAGS) -c $(INCDIRS) -o $@ $<

$(ODIR)/%.obj : %.c $(HDRS) $(ODIR)
	$(CC) $(CFLAGS) -c $(INCDIRS) -o $@ $<

$(ODIR):
	mkdir $@
