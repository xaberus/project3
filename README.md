project3
========

These are the sources for projekt number 3 for computational physics course I am participating.
The purpose of this project is to appoximate solutions to the time dependant Schrödinger-equation
and retreive a spectrum of stationary energy eigenvalues by calculating the autocorrlation
    (psi|U(t)|psi)
using an approximated unitary time propagator U.


To start the simulation a configuration file is needed, which is parsed by the means of Lua 5.1, i.e.

    config = {
      bins = bins; dt = 5.7/1000;
      range = {-dx*bins/2,dx*bins/2};
      steps = 1; runs = 16384;
      potential = function(x)
        local k0, k2, k3, k4 = -132.7074997, 7, .5, 1
        return  k0-k2*x^2+k3*x^3+k4*x^4
      end;
      psi = function(x)
        local a, s = 1.9, 0.87
        return {exp(-(x-a)^2/(2*s^2)),exp(-(x+a)^2/(2*s^2))}
      end;
      enrgrange = {-139, -137.8, 6, 1.5};
      output = {
        dir = "./badness";
        apsi = "apsi.dat";
        pot = "pot.dat";
        corr = "corr.dat";
        dftcorr = "dftcorr.dat";
        theoenrg = "theoenrg.dat";
        spectrum = "spectrum.dat";
      };
    }


The programm consists of two main parts.
* simulate
* evaluate

simulate
--------

    simulate <config.lua> <results.dat>
    
simulate takes a configuration file and runs the simulation according to the specified parameters.
The results are written to the results file afterwards.

evaluate
--------

    evaluate <config.lua> <results.dat>

evaluate takes both the config and the results file and evaluates the results so they can be viewed
in a human readable form. evaluate also conducts a peak search and an approximation by the Numerov method,
to have comparable values.

licence
=======

All the source code and documentation written by me is aviable under the MIT licence, but as this project 
uses fftw3 for it's calculations the overall licence must be GPLv2. (As long as you replace the fftw
parts, this project my be regarded as MIT-licenced.)

complex.lua is taken verbatim from http://lua-users.org/wiki/ComplexNumbers and is also MIT licenced as
is Lua itself.