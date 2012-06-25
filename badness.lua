local math = require("math")
local complex = require("complex")
local exp = math.exp

--local dx = 0.825/10
--local bins = 512
local bins = 4096
local dx = 16/bins

config = {
  bins = bins;
  dt = 5.7/1000;
  range = {-dx*bins/2,dx*bins/2};
  steps = 1; runs = 16384;
  --
  potential = function(x)
    --local k0, k2, k3, k4 = -100, 40, 2, 3
    local k0, k2, k3, k4 = -132.7074997, 7, .5, 1
    return  k0-k2*x^2+k3*x^3+k4*x^4
  end;
  psi = function(x)
    local a, s = 1.9, 0.87
    return {exp(-(x-a)^2/(2*s^2)),exp(-(x+a)^2/(2*s^2))}
  end;
  energy = function(k)
    return 1/0
  end;
  enrgrange = {-150, 0, 6};
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
