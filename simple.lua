local math = require("math")
local complex = require("complex")
local exp = math.exp

config = {
  bins = 4096*4;
  dt = 0.0001;
  range = {-8,8};
  steps = 10; runs = 100000;
  --steps = 10; runs = 100000;
  vstep = 1; vframes = 1;
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
    local x = k + 1
    return
      -2.5805486322214
      + 9.1857765367527 * x
      + 0.08881127950008 * x^2
      - 0.000260375459505 * x^3
  end;
  enrgrange = {0, 205, 6};
  output = {
    dir = "./simple";
    apsi = "apsi.dat";
    pot = "pot.dat";
    corr = "corr.dat";
    dftcorr = "dftcorr.dat";
    theoenrg = "theoenrg.dat";
    spectrum = "spectrum.dat";
  };
}