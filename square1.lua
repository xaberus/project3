local math = require("math")
local complex = require("complex")
local L = 10
local exp = math.exp
local cexp = complex.exp
local pi = math.pi
local sqrt = math.sqrt
local sin = math.sin

config = {
  bins = 4096;  dt = 0.0001;
  range = {-5,15};
  steps = 10; runs = 50000;
  vstep = 100; vframes = 200;
  --
  potential = function(x)
    if 0 < x and x < L then return 0 else return 10000 end
  end;
  psi = function(x)
    if 0 < x and x < L then
      return {
        -- + sqrt(2/L) * sin(2*pi/L * x)
        sqrt(2/L) * sin(1*pi/L * x),
        0
      }
    else
      return  {0,0}
    end
  end;
  energy = function(r)
    local k = r + 1
    return pi^2 * k^2/L^2
  end;
  enrgrange = {0, 1, 8, 2};
  output = {
    dir = "./square1";
    apsi = "apsi.dat";
    pot = "pot.dat";
    corr = "corr.dat";
    dftcorr = "dftcorr.dat";
    spectrum = "spectrum.dat";
  };
}