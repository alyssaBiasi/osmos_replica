UNAME := $(shell uname)
G = g++

CFLAGS = -ansi -Wall -pedantic -c -g
LDFLAGS = `sdl-config --libs` -lglut -lGLU -lGLEW -lGL -lm
EXE = ass3
OBJECTS = assignment3.o
 
# OS X
ifeq "$(UNAME)" "Darwin"
	G = g++
	LDFLAGS = -framework Carbon -framework OpenGL \
			-framework GLUT -lz
endif

#default target
all: $(EXE)

#executable
$(EXE) : $(OBJECTS)
	$(G) -o $@ $(LDFLAGS) $(OBJECTS)

#general object file rule
%.o : %.c
	$(G) -o $@ $(CFLAGS) libSOIL.a $<

#additional dependencies
particles_2D.o : utils.h

#clean (-f stops error if files don't exist)
clean:
	rm -f $(EXE) $(OBJECTS)

