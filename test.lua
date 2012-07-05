local math = require("math")
local complex = require("complex")
local omega = 5;
local exp = math.exp
local sqrt = math.sqrt
local pi = math.pi

config = {
  bins = 4096; dt = 1 / (1000 * pi);
  range = {-15,15};
  steps = 10; runs = 50000;
  vstep = 100; vframes = 200;
  --
  potential = function(x)
    return omega^2/4 * x ^ 2;
  end;
  psi = function(x)
    local a = .5
    local aa = a * a
    local x0 = 0
    local k0 = 10
    local xx = (x-x0) * (x-x0)
    return  math.exp(-xx/(aa)) * complex.exp({0, k0*(x)})
  end;
  energy = function(k)
    return omega * (.5 + k)
  end;
  enrgrange = {0, 200, 6, 1};
  output = {
    dir = "./test";
    apsi = "apsi.dat";
    pot = "pot.dat";
    corr = "corr.dat";
    dftcorr = "dftcorr.dat";
    spectrum = "spectrum.dat";
  };
}