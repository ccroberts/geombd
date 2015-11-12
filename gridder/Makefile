CC = icpc
CFLAGS = -O3 -std=c++0x -I/usr/include/x86_64-linux-gnu/c++/4.8/ -openmp -I../
#CFLAGS += -DDEBUG
#CFLAGS += -DTOUCHING
LDFLAGS = -lcilkrts -openmp
		
O_gridder = Main.o 
X_gridder = gridder
		
O_griddes = ES.o 
X_griddes = griddes

O_spm2dx = spm2dx.o 
X_spm2dx = spm2dx

O_bpm2dx = bpm2dx.o 
X_bpm2dx = bpm2dx

O_griddex = Exclusion.o
X_griddex = griddex

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
