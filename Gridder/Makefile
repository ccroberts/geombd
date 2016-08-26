CC = icpc
CFLAGS = -O3 -std=c++0x -I/usr/include/x86_64-linux-gnu/c++/4.8/ -openmp -I../GeomBD
#CFLAGS += -DDEBUG
#CFLAGS += -DTOUCHING
LDFLAGS = -lcilkrts -openmp
		
O_gLJ = Gridder-LJ.o 
X_gLJ = ../bin/Gridder-LJ
		
O_gEX = Gridder-EX.o 
X_gEX = ../bin/Gridder-EX
		
O_gES = Gridder-ES.o 
X_gES = ../bin/Gridder-ES
		
O_gD = Gridder-D.o 
X_gD = ../bin/Gridder-D
		
O_gBPM2DX = bpm2dx.o
X_gBPM2DX = ../bin/bpm2dx
		
O_gBXM2DX = bxm2dx.o
X_gBXM2DX = ../bin/bxm2dx

all : gLJ gEX gES gD gBPM2DX gBXM2DX

gLJ: $(O_gLJ) $(C_gLJ)
	$(CC) $(CFLAGS) $(O_gLJ) $(C_gLJ) -o $(X_gLJ) $(LDFLAGS)

gEX: $(O_gEX) $(C_gEX)
	$(CC) $(CFLAGS) $(O_gEX) $(C_gEX) -o $(X_gEX) $(LDFLAGS)

gES: $(O_gES) $(C_gES)
	$(CC) $(CFLAGS) $(O_gES) $(C_gES) -o $(X_gES) $(LDFLAGS)

gD: $(O_gD) $(C_gD)
	$(CC) $(CFLAGS) $(O_gD) $(C_gD) -o $(X_gD) $(LDFLAGS)

gBPM2DX: $(O_gBPM2DX) $(C_gBPM2DX)
	$(CC) $(CFLAGS) $(O_gBPM2DX) $(C_gBPM2DX) -o $(X_gBPM2DX) $(LDFLAGS)

gBXM2DX: $(O_gBXM2DX) $(C_gBXM2DX)
	$(CC) $(CFLAGS) $(O_gBXM2DX) $(C_gBXM2DX) -o $(X_gBXM2DX) $(LDFLAGS)

%.o : %.cc
	@echo " - Compiling $<..." 
	$(CC) $(CFLAGS) -c $<

%.o : %.cu
	@echo " - Compiling $<..." 
	$(NVCC) -c $<

clean:
	@rm *.o $(X_gLJ) $(X_gES) $(X_gD) $(X_gBPM2DX) $(X_gBXM2DX) $(X_gEX)
