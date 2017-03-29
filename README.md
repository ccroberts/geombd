# geombd 3.0beta


GeomBD is a rigid body Brownian dynamics software for determining interenzyme intermediate transfer rates and substrate association rates in biomolecular complexes. GeomBD2 supports the <a href="http://dx.doi.org/10.1063/1.446900">NAM method</a> for determining second-order association rate constants. In addition, GeomBD2 extends this scheme to quantify intermediate substrate transfer rates in multi-enzyme complexes. Thus, simulations in GeomBD2 can proceed through two different configurations: ligand-receptor association, and substrate transfer. For both types of simulations, many substrate replicates diffuse in parallel among the spatially fixed, static descriptions of the receptor system. 
<br>

The following features are included in GeomBD2:

* Simple rigid body ligand PQR-formatted definitions composed of charged spherical particles.
* Volume exclusion grids for receptor representation.
* Precalculated electrostatic grid support for including electrostatic forces in simulations via the screened Coulomb potentia.
* Precalculated Lennard-Jones potential grid support for short-range forces and receptor representation.
* Linearly variable timestep.
* Reaction criteria defined by distance threshold between ligand atom(s) and binding site coordinate(s).
* Traditional NAM method for association rate constant determination.
* Inter-receptor transfer simulations using a split scheme determining "direct" and "indirect" contributions to intermediate transfer rate
* Automated termination of simulation sets via convergence criteria.
* Intel Cilk+ threaded parallel simulation.

### Requirements
* Linux or Unix-like operating system
* Intel C++ compiler (or GCC with support for CilkPlus)
* Python 2.6+ for supporting analysis scripts

### Input File
* __temperature__ [Kelvin]:<br> Specify the temperature of the system<br><br>
* __timestep__ [Smallest Timestep (ps)] [Largest Timestep (ps)] [Radial scale start (A)] [Radial scale end (A)]:<br> Specify the linearly scaling timestep. The first two values specificy the minimum and maximum timestep values, while the last two values specify the radial distances from which the timestep is scaled from the minimum to the maximum value.<br><br>
* __order__ [integer]:<br> Specify the order of finite difference force approximation for potential grids. (Default=2)<br><br>
* __threads__ [# of threads]:<br> Specify the number of the parallel threads (Cilk workers) to use<br><br>
* __writetraj__ [# of steps]:<br> Specify the frequency of trajectory writes in number of steps per write.<br><br>
* __writelog__ [# of steps]:<br> Specify the frequency of reporting rate constant data to log file<br><br>
* __logbinders__:<br> If this keyword appears in the input file, individual binding events will be logged.<br><br>
* __logexiters__:<br> If this keyword appears in the input file, unsuccessful simulation replicates that leave the simulation space will be logged.<br><br>
* __convergence__ [Convergence criteria]:<br> Specify the criteria with which the convergence of the Beta value (bound ligand fraction) is determined. The convergence test is based on the standard error of the Beta value. As an example, a __convergence__ value of 0.0001 will ensure that the standard error of Beta is <= 0.0001.<br><br>
* __convwindow__ [# of completed simulations]:<br> Specify the number of the most recent completed replicate simulations to consider when calculating the standard error for the convergence test.<br><br>
* __receptor__ [PQR filename]:<br> Specify the receptor system for the simulation. This file is used to determine the center of the simulation space. Actual representations of the receptor system must be specified by grid files.<br><br>
* __grid__ [Grid Type] [BPM/BXM grid filename]:<br> Specify a grid-based representation of a receptor.<br>*Grid Type* can be:<br>
    * __ex__ for an exclusion grid (BXM)<br>
    * __es__ for an electrostatic potential grid (BPM)<br>
    * __d__ for a desolvation potential grid (BPM)<br>
    * a ligand atom type as it appears in the ligand PQR file (ex:<br> C, O, H, N) for a Lennard-Jones potential grid (BPM).<br><br>
* __session__ [direct / indirect / nam]:<br> Begin definition of simulation session block. __direct__ specifies a first-order intermediate transfer simulation, while __indirect__ and __nam__ are synonymous for a second-order association rate determination.<br><br>
    * __ligand__ [PQR filename]:<br> Specifies the structure file for a ligand to be diffused in the current simulation _session_ definition. This can be a single structure or multiple ligand conformations separated by END statements.<br><br>
    * __b__ [radius, Angstroms]:<br> Specifies the b-radius for the system. For NAM/indirect simulations, this is the starting radius of all ligands. For __direct__ transfer simulation types, this is the termination radius.<br><br>
    * __q__ [radius, Angstroms]:<br> Specifies the q-radius (termination radius) for the __indirect/NAM__ simulations.<br><br>
    * __from__ [X] [Y] [Z]:<br> Specifies the starting coordinate for ligands in __direct__ simulations.<br><br>
    * __bind__ [X] [Y] [Z] [Radius]:<br> Specifies the binding site cartesian coordinate and binding site radius, in Angstroms.<br><br>

Example intermediate transfer input file:
<div style="border: 1px solid #777; background-color: #f0f0f0; padding: 10px;">
  <code>
  temperature 298.15<br>
  timestep 0.050 0.500 75 300<br>
  threads 20<br>
  writetraj 10000<br>
  writelog 10000<br>
  convergence 0.0001<br>
  convwindow 1000<br>
  <br>
  receptor H95.pqr<br>
  grid ex H95-ex.bxm<br>
  <br>
  session direct<br>
  &nbsp;&nbsp;ligand Intermediate.pqr 1600<br>
  &nbsp;&nbsp;b 150<br>
  &nbsp;&nbsp;from 57.092442 31.592548 65.040970<br>
  &nbsp;&nbsp;bind 39.2522 43.8311 43.4662 7 5.0<br>
  <br>
  session indirect<br>
  &nbsp;&nbsp;ligand Intermediate.pqr 1600<br>
  &nbsp;&nbsp;b 150<br>
  &nbsp;&nbsp;q 700<br>
  &nbsp;&nbsp;bind 39.2522 43.8311 43.4662 7 5.0<br>
  </code>
</div>



GeomBD can be run using the following example command:
<div style="border: 1px solid #777; background-color: #f0f0f0; padding: 10px;">
  <code>
    GBD -i INPUT_FILE -l LOG_FILE -o TRAJ_FILE.pqr
  </code>
</div>


<a href="https://www.lucidchart.com/invitations/accept/994b3d12-de47-4c70-a02c-8b096711f0e7">Basic source code diagram</a>
