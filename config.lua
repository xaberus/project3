local math = require("math")
local complex = require("complex")

local omega = 40;
local A = math.sqrt(4 * omega)

config = {
  --bins = 2 * 4096;
  bins = 128;
  dt = 0.0001;
  range = {-25,25};
  steps = 10;
  runs = 10000;
  --
  vstep = 100;
  vframes = 200;
  --
  potential = function(x)
    return omega * x ^ 2;
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
    return A * (.5 + k)
  end;
  output = {
    dir = "./harmosca";
    corr = "corr.dat";
    dftcorr = "dftcorr.dat";
  };
}