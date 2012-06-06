local math = require("math")
local complex = require("complex")
local omega = 40;
local L = 10
local exp = math.exp
local sin = math.sin
local sqrt = math.sqrt
local cexp = complex.exp
local pi = math.pi

config = {
  bins = 4096*2;  dt = 0.0001;
  range = {-10,20};
  steps = 10; runs = 100000;
  vstep = 100; vframes = 200;
  --
  potential = function(x)
    if 0 < x and x < L then return 0 else return 10000 end
  end;
  psi = function(x)
    return exp(-(x-L/2)^2/(2)) * cexp({0, 2 * x})
  end;
  energy = function(k)
    return pi^2 * k^2/L^2
  end;
  enrgrange = {0, 20, 4};
  output = {
    dir = "./square";
    apsi = "apsi.dat";
    pot = "pot.dat";
    corr = "corr.dat";
    dftcorr = "dftcorr.dat";
    theoenrg = "theoenrg.dat";
    spectrum = "spectrum.dat";
  };
}