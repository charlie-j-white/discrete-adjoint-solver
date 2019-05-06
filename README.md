# Discrete Adjoint Solver

By Charlie Anderson for the UoB MEng Aerospace Engineering Final Year Project module






## Getting started

If you have not downloaded this repository before, you will need to use git to clone a new one. Just copy and paste this into the command line:
```
$ git clone https://github.com/charlie-j-white/discrete-adjoint-solver.git
```
Then perform the following instructions to get a working executable:






### Dependencies

The mandatory dependencies are as follows -

- `gfortran`
- `gcc`
- LAPACK

For Arch-based distributions LAPACK can be downloaded with
```
$ sudo pacman -S lapack
```
or for Ubuntu
```
$ sudo apt-get install liblapack-dev
```

##### Optional dependencies:

The optional dependencies are -

- `python`
- Matplotlib

in the case where you do not have Tecplot installed and need some way of visualising the program outputs. The Python scripts in `DEV/` use Matplotlib to read the `.plt` files, which can be installed with
```
$ sudo pip install matplotlib
```





### Build instructions

This program can easily be compiled on an up to date Linux system. Just follow the steps
```
$ cd discrete-adjoint-solver/
$ make
$ ./main
```
and this will create and run the executable. Feel free to edit the Makefile if needed.













## Program use

For an in-depth view into the theory and motivation behind the code, refer to the report, but here the process of using the code will be explained.


### Configuration

Because automatic differentiation does not work well with `allocate()` dynamic memory statements this program will have to be recompiled each time a parameter is changed. Fortunately this does not take very long as all of the source code files are compiled individually before being linked.

All the changes to be made should be done editing the
```
$ edit SRC/wrapper.f
``` 
function. While making changes to the other functions is permitted, this is for development purposes only, and of course your fault if you then break it. Unless otherwise stated, all configuration instructions refer to this file.




##### Runtype

On line `36` the runtype option is available, with the options listed below. This specifies the task performed by the program:

```
     runtype = 0: Runs the primal CFD code.

     runtype = 1: Sensitivity via finite differences

     runtype = 2: Sensitivity via the adjoint method

     runtype = 3: Runs the adjoint method next to finite differences for comparison
```

##### Mesh size

On lines `77` and `78` are the variables controlling the number of x and y cells for the structured mesh. As the current linear solver deals with the badly conditioned flux Jacobian the adjoint results get worse with increasing mesh size. For values of `nx=20` and `ny=3` close comparisons with finite differences are output, though of course this would handle bigger meshes if the linear solver was replaced.


##### Design variables

In order to fully demonstrate the adjoint method the number of design variables was designed to be changed easily. As the shape is initialised in the `wrapper` function this must be edited from the variable declarations. On line `8` the variable `na` is the number of design variables, but **this must be changed at the same time** as the array size on line `11` in order to have the right size.



##### Other parameters

Most flow instructions are contained in the `params` variable, defined from line `87` onwards. This is fairly self-explanatory, but change these at your own risk.


The shape of the mesh can be changed as well. On line `137` the mesh is initialised in a cos^2 wave shape, but this can be changed to anything. Bear in mind that the start and end point of the duct height are fixed, so this cannot change.


##### Flux Jacobian visualisation

The entire program is written in FORTRAN, with the exception of the external library `LODEPNG` used to create an image of the flux Jacobian. By default, this outputs an image unless the file
```
$ edit SRC/adjoint.f
```
is edited and the line
```
cerr = colour(nw,fluxjac)
```
is commented out.



 

### Output files

Different runtypes give different output files. A few will be listed here:

- `solution.plt` is the full two-dimensional flow solution which should be used for input into Tecplot

- `convergence.plt` shows the RMS convergence against iteration, used by the `DEV/convdisp.py` script

- `sensitivities.plt` lists the sensitivities calculated by both finite differences and the adjoint method. Best used with `runtype=3` but works fine with `=2` or `=3`. Remember this works best with small mesh sizes due to matrix conditioning. Used by `DEV/sensdisp.py`

- `memory.plt` shows the RSS memory in kB of the program over time as defined by the status file in `/proc/[pid]/status`. This program is created by the auxiliary `./memory` executable which can be compiled in the `DEV/` folder, and is read by the `memdisp.py` script

- `transect.plt` is the solution taken down the centre of the duct, used of quasi-1D flow visualisation by `DEV/display.py`

- `flux-jacobian.png` a PNG image of the flux Jacobian. All non-zero elements are shown in black



### Developer files

Developer files, stored in the `DEV/` directory have been created. As well as for debugging purposes, they are designed to process the data rather than refer to third-party software. Note that all of these are designed to be **run from the top directory**, meaning they will need to be moved there first. Try something like
```
$ mv display.py ../
```


The Python scripts (`DEV/*.py`) are used to read the solution files and display the data using Matplotlib. The scripts needed to read each solution file are shown above. They can be run by using
```
$ python display.py
```
Note that no further arguments are needed.



The shell scripts `*.sh` are used for various things, mainly for redundant debugging and development processes.


The C program used to create an image of the flux Jacobian is as previously described.



The final developer program is the one designed to measure memory usage. This is a program written in FORTRAN. To build it:
```
$ cd DEV/
$ make
$ mv memory ../
$ cd ../
```
At this stage, you should be in the top directory. The main program you want to measure the memory usage of should be already built, so binary executables `main` and `memory` should be there. To use the memory program, type
```
./memory
```
in one terminal. This should start searching for the main program. Then in a different terminal, type
```
./main
```
which will run the adjoint solver. At this stage the memory program should recognise that the main program has started running, and start measuring memory usage. When they have both concluded, output `memory.plt` will be outputted.





### Troubleshooting

When using the Python scripts, if you get the `wrong number of columns` error this means the numbers are too large in the relevant solution file, and may need to be separated manually. Or the code needs to be re-run with smaller programs. Either way, this is something I should probably fix at some point.












