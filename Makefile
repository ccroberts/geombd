
all:
	cd GeomBD; make
	cd Gridder; make

clean:
	cd GeomBD; make clean
	cd Gridder; make clean
