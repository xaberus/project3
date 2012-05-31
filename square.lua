local math = require("math")
local complex = require("complex")

local omega = 40;
local L = 10

config = {
  bins = 4096*2;
  dt = 0.0001;
  range = {-25,25};
  steps = 10;
  --runs = 10000;
  runs = 50000;
  --
  --vstep = 100;
  --vframes = 200;
  --
  potential = function(x)
    if math.abs(x) < L/2 then return 0 else return 10000 end
  end;
  psi = function(x)
    return  math.exp(-x^2/(2)) * complex.exp({0, 4 * x})
  end;
  energy = function(k)
    return A * (.5 + k)
  end;
  output = {
    dir = "./square";
    apsi = "apsi.dat";
    pot = "pot.dat";
    corr = "corr.dat";
    dftcorr = "dftcorr.dat";
  };
}