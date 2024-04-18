# Makefile

TOP_SRCDIR = /home/demulab/src/ode-0.16.1

CXX = g++
CFLAGS = -g

LIB = -lode -ldrawstuff -lGL -lGLU -lX11 -lrt -lm -lpthread -lstdc++
LIBDIR = -L/usr/lib/x86_64-linux-gnu -L$(TOP_SRCDIR)/drawstuff/src/.libs  -L$(TOP_SRCDIR)/ode/src/.libs -L/usr/lib
INCDIR = -I$(TOP_SRCDIR)/include -I$(TOP_SRCDIR)/ode/src -I/usr/include


OBJS =xythetaLearning action_read stickness2

all: $(OBJS)

$(OBJS): %: %.cpp texturepath.h
	$(CXX) $(CFLAGS) -o $@ $< $(INCDIR) $(LIBDIR) $(LIB)

clean:
	rm  $(OBJS)  *~  *.*~


.PHONY: all clean
