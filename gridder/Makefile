CC = icpc
CFLAGS = -O3 -std=c++0x -I/usr/include/x86_64-linux-gnu/c++/4.8/ -openmp -I../
#CFLAGS += -DDEBUG
#CFLAGS += -DTOUCHING
LDFLAGS = -lcilkrts -openmp
		
O_gridder = Main.o 
X_gridder = gridder

O_bpm2dx = bpm2dx.o 
X_bpm2dx = bpm2dx

all : gridder bpm2dx

gridder: $(O_gridder) $(C_gridder)
	$(CC) $(CFLAGS) $(O_gridder) $(C_gridder) -o $(X_gridder) $(LDFLAGS)

bpm2dx: $(O_bpm2dx) $(C_bpm2dx)
	$(CC) $(CFLAGS) $(O_bpm2dx) $(C_bpm2dx) -o $(X_bpm2dx) $(LDFLAGS)

%.o : %.cc
	@echo " - Compiling $<..." 
	$(CC) $(CFLAGS) -c $<

%.o : %.cu
	@echo " - Compiling $<..." 
	$(NVCC) -c $<

clean:
	@rm *.o $(X_gridder)
