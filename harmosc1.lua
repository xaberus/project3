local math = require("math")
local complex = require("complex")
local omega = 6;
local exp = math.exp
local sqrt = math.sqrt

config = {
  bins = 4096*2; dt = 0.0001;
  range = {-20,20};
  steps = 10; runs = 50000;
  vstep = 100; vframes = 200;
  --
  potential = function(x)
    return omega^2/4 * x ^ 2;
  end;
  psi = function(x)
    --return  {exp(-.5*omega/4/sqrt(2.5)*x^2),0}
    --return  {exp(-omega/4*x^2),0}
    return {
      exp(-omega/4*x^2) + (x*sqrt(omega))/exp(x^2*omega/4.),
      0
    }
  end;
  energy = function(k)
    return omega * (.5 + k)
  end;
  enrgrange = {0, 2*omega, 6};
  output = {
    dir = "./harmosc1";
    apsi = "apsi.dat";
    pot = "pot.dat";
    corr = "corr.dat";
    dftcorr = "dftcorr.dat";
    theoenrg = "theoenrg.dat";
    spectrum = "spectrum.dat";
  };
}