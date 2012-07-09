local math = require("math")
local complex = require("complex")
local exp = math.exp

config = {
  cmp = true;
  bins = 4096;
  dt = 0.0001;
  range = {-8,8};
  steps = 10; runs = 150000;
  vstep = 1; vframes = 1;
  --
  potential = function(x)
    local k0, k2, k3, k4 = -132.7074997, 7, .5, 1
    return  k0-k2*x^2+k3*x^3+k4*x^4
  end;
  psi = function(x)
    local a, s = 1.9, 0.87
    return {exp(-(x-a)^2/(2*s^2)),exp(-(x+a)^2/(2*s^2))}
  end;
  enrgrange = {-146, 0, 4, 1};
  output = {
    dir = "./data/cmp";
    apsi = "apsi.dat";
    pot = "pot.dat";
    corr = "corr.dat";
    dftcorr = "dftcorr.dat";
    spectrum = "spectrum.dat";
    numen = "numen.dat";
    splen = "splen.dat";
    aken = "aken.dat";
  };
}
