local math = require("math")
local complex = require("complex")
local L = 10
local exp = math.exp
local cexp = complex.exp
local pi = math.pi

config = {
  bins = 4096;  dt = 0.0001;
  range = {-5,15};
  steps = 10; runs = 50000;
  vstep = 100; vframes = 200;
  --
  potential = function(x)
    if 0 < x and x < L then return 0 else return 1000 end
  end;
  psi = function(x)
    return exp(-(x-L/2)^2/(2)) * cexp({0, 2 * x})
  end;
  energy = function(r)
    local k = r + 1
    return pi^2 * k^2/L^2
  end;
  enrgrange = {0, 20, 4, 2};
  output = {
    dir = "./data/square";
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