local math = require("math")
local complex = require("complex")

local exp = math.exp

config = {
  bins = 4096;
  dt = 0.0001;
  range = {-7,7};
  steps = 10;
  runs = 20000;
  --
  --vstep = 100;
  --vframes = 200;
  --
  potential = function(x)
    local k0, k2, k3, k4 = -100, 20, .5, 1
    return  k0-k2*x^2+k3*x^3+k4*x^4
  end;
  psi = function(x)
    local a, s = 1.9, 0.87
    return exp(-(x-a)^2/(2*s^2)) + exp(-(x+a)^2/(2*s^2)) * complex.exp({0, 4 * x})
  end;
  energy = function(k)
    return k*2000
  end;
  enrgrange = {0, 205};
  output = {
    dir = "./simple";
    apsi = "apsi.dat";
    pot = "pot.dat";
    corr = "corr.dat";
    dftcorr = "dftcorr.dat";
    theoenrg = "theoenrg.dat";
  };
}
