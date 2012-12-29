Cosmocal
=========

Cosmocal is a c++ library package that computes different cosmological quantities

It calculates cosmic time, ang. diameter distance, growth factor and its log derivative and few other cosmological quantities.
input- cosmological parameters and redshift

Cosmocal is released under the MIT software liscense (see LISCENSE).

Prerequisites
=============

Gnu Scientific Library (GSL)
A c++ compiler

Descriptions
============

The package contain 3 header files- cosmo.h, growth.h and driver.h, 1 input parameter file- input.par, one makefile-Makefile and one driver routine- cosmo_driver.cpp

1. cosmo.h- contains the class "cosmo" and the associated functions that calculates most of the cosmogical quantitites e.g., angular diamter distance and cosmic time.
2. growth.h- contain the class "growth" and calulates the growth factor and its log derivative. (N.B.- does not assume any approximate form for the growth factor and solve the ODE )
3. driver.h- defines the input parameters
4. Makefile- makefile to generate the executable
5. input.par- file that contain the values of the input parameters.
6. cosmo_driver.cpp- the main driver routine that calls the class cosmo and growth and computes the cosmological quantities

 

Installation and running
========================

1. Download the package and place it anywhere you like
2. Open the Makefile, edit CFLAGS and CLIB to make sure you have the correct path to GSL insalled in your machine.
3. do make clean, then make cosmo, this creates the executable
4. run : ./cosmo_driver.exe < input.par
