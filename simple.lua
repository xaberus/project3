local math = require("math")
local complex = require("complex")

local exp = math.exp

--[=====================[
tab = Import["/tmp/projekt/simple/dftcorr.dat", "Table"];
use = 1000;
xytab = {#[[1]], #[[4]]} & /@ tab[[1 ;; use]];
xtab = #[[1]] & /@ tab[[1 ;; use]];
ytab = #[[4]] & /@ tab[[1 ;; use]];
nytab = (ytab - Min[ytab])/(Max[ytab] - Min[ytab]);
maxs = Floor[Median@#] & /@
   Split[Most[
      ArrayRules[SparseArray[If[# > .02, 1, 0] & /@ nytab]]] /.
     HoldPattern[{a_} -> _] :> a, Abs[#2 - #1] < 7 &];
unga = Transpose[{maxs, nytab[[maxs]]}];
xmpos = xtab[[maxs]];
preg = Table[{n, xmpos[[n]]}, {n, 1, Length[xmpos]}];
reg = Fit[preg, {x^4, x^2, x, 1}, x]
card = Table[{Block[{x = k}, reg], 0.000025}, {k, 1, Length[preg]}];
Show[ListLinePlot[xytab, PlotRange -> {Automatic, {0, 0.00003}}],
 ListPlot[card]]
Show[ListPlot[{preg,
   Table[Block[{x = k}, reg], {k, 1, Length[preg]}]}],
 Plot[reg, {x, 0, 25}]]
]=====================]

config = {
  bins = 4096;
  dt = 0.0001;
  range = {-8,8};
  steps = 10;
  runs = 50000;
  --
  vstep = 100;
  vframes = 200;
  --
  potential = function(x)
    local k0, k2, k3, k4 = -100, 40, .5, 1
    return  k0-k2*x^2+k3*x^3+k4*x^4
  end;
  psi = function(x)
    local a, s = 1.9, 0.87
    return complex({0,-exp(-(x-a)^2/(2*s^2))}) + exp(-(x+a)^2/(2*s^2)) * complex.exp({0, 4 * x})
  end;
  energy = function(k)
    local x = k + 1
    return -5.2527429168021 + 6.9743356909693 * x + 0.1675516081915 * x^2 - 0.004188790204786 * x^4
    --return -4.2818872178543 + 6.8630652952975 * x + 0.06073125879362 * x^2 - 0.000410250052803 * x^3
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
