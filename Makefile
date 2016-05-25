
all:
	cd GeomBD; make
	cd Gridder; make
	cd ProbDX; make

clean:
	cd GeomBD; make clean
	cd Gridder; make clean
	cd ProbDX; make clean
