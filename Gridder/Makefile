CC = icpc
CFLAGS = -O3 -std=c++0x -I/usr/include/x86_64-linux-gnu/c++/4.8/ -openmp -I../GeomBD
#CFLAGS += -DDEBUG
#CFLAGS += -DTOUCHING
LDFLAGS = -lcilkrts -openmp
		
O_gridder = Gridder-LJ.o 
X_gridder = ../bin/Gridder-LJ
		
O_griddes = Gridder-ES.o 
X_griddes = ../bin/Gridder-ES

O_spm2dx = spm2dx.o 
X_spm2dx = ../bin/Gridder-spm2dx

O_bpm2dx = bpm2dx.o 
X_bpm2dx = ../bin/Gridder-bpm2dx

O_griddex = Gridder-EX.o
X_griddex = ../bin/Gridder-EX

all : gridder griddes griddex bpm2dx spm2dx

gridder: $(O_gridder) $(C_gridder)
	$(CC) $(CFLAGS) $(O_gridder) $(C_gridder) -o $(X_gridder) $(LDFLAGS)

griddes: $(O_griddes) $(C_griddes)
	$(CC) $(CFLAGS) $(O_griddes) $(C_griddes) -o $(X_griddes) $(LDFLAGS)

griddex: $(O_griddex) $(C_griddex)
	$(CC) $(CFLAGS) $(O_griddex) $(C_griddex) -o $(X_griddex) $(LDFLAGS)

bpm2dx: $(O_bpm2dx) $(C_bpm2dx)
	$(CC) $(CFLAGS) $(O_bpm2dx) $(C_bpm2dx) -o $(X_bpm2dx) $(LDFLAGS)

spm2dx: $(O_spm2dx) $(C_spm2dx)
	$(CC) $(CFLAGS) $(O_spm2dx) $(C_spm2dx) -o $(X_spm2dx) $(LDFLAGS)

%.o : %.cc
	@echo " - Compiling $<..." 
	$(CC) $(CFLAGS) -c $<

%.o : %.cu
	@echo " - Compiling $<..." 
	$(NVCC) -c $<

clean:
	@rm *.o $(X_gridder) $(X_griddes) $(X_bpm2dx) $(X_griddex)
