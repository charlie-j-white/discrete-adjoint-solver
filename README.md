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
and this will create an executable. Feel free to edit the Makefile if needed.




## Program use


How to use the program



