
local dump = require "dump"

local tol = .3

local abs = math.abs

local args = {...}
local sep = "\%s+"

function split(self, inSplitPattern, outResults )
  if not outResults then
    outResults = { }
  end
  local theStart = 1
  local theSplitStart, theSplitEnd = string.find( self, inSplitPattern,
theStart )
  while theSplitStart do
    table.insert( outResults, string.sub( self, theStart, theSplitStart-1 ) )
    theStart = theSplitEnd + 1
    theSplitStart, theSplitEnd = string.find( self, inSplitPattern, theStart )
  end
  table.insert( outResults, string.sub( self, theStart ) )
  return outResults
end


local function read_table(file)
  io.stderr:write(string.format("loading '%s'\n", file))
  local fp = assert(io.open(file, "r"), string.format("could not open %s", file))
  local reader = fp:lines()
  local ftab = {}
  while true do
    local line = reader()
    if line then
      if line:match("^[^#]") then
        local res = split(line, sep)
        ftab.cols = #res
        ftab[#ftab+1] = res
      end
    else
      break
    end
  end
  return ftab
end

assert(#args == 3, "usage: tabgen config.lua spectrum.dat numen.dat")
loadfile(args[1])() -- config
assert(config)
local pre_spectab = read_table(args[2])
local pre_numtab = read_table(args[3])

local vlst = {}

local function nconcat(prow)
  local row = {}
  for k, v in ipairs(prow) do
    local t = type(v)
    if t == "number" then
      for k = 1, v, 1 do
        row[#row+1] = ""
      end
    elseif t == "table" then
      for j, w in ipairs(v) do
        row[#row+1] = w
      end
    end
  end
  return row
end

local function tonumbers(row)
  local nrow = {}
  for k, v in ipairs(row) do
    nrow[#nrow+1] = tonumber(v)
  end
  return nrow
end

local enfn = config.energy
local entab
if enfn then
  local Emax = assert(config.enrgrange[2])
  entab = {}
  local k = 0
  local E = 0
  repeat
    E = enfn(k)
    if E < Emax then
      entab[#entab + 1] = {E}
      k = k + 1
    end
  until E > Emax
  entab = entab
end

local samp = {pre_spectab, pre_numtab, entab}
local indices = {1,1,1}

while true do
  local irow = {}
  local empty = true
  for k, s in ipairs(samp) do
    local ss = s[indices[k]]
    if ss then
      irow[k] = tonumbers(ss)
      empty = false
    else
      irow[k] = nil
    end
  end

  if empty then
    break
  else
    local min
    for k, ss in pairs(irow) do
      if min then
        if ss[1] < min then
          min = ss[1]
        end
      else
        min = ss[1]
      end
    end
    for k, ss in pairs(irow) do
      local d = (min - ss[1])
      if abs(d) < tol then
        indices[k] = indices[k] + 1
      else
        irow[k] = nil
      end
    end

    for k=1, #samp, 1 do
      if not irow[k] then
        irow[k] = samp[k].cols or 1
      end
    end

    vlst[#vlst+1] = nconcat(irow)
  end
end

--dump[1](vlst)

--[[
local i, j = 1, 1
while true do
  local er, nr = pre_spectab[i], pre_numtab[j]
  if er and nr then
    er = tonumbers(er)
    nr = tonumbers(nr)
    local d = (er[1] - nr[1])
    if abs(d) < tol then
      vlst[#vlst+1] = nconcat({er, nr})
      i = i + 1
      j = j + 1
    elseif d > 0 then
      vlst[#vlst+1] = nconcat({pre_spectab.cols, nr})
      j = j + 1
    elseif d < 0 then
      vlst[#vlst+1] = nconcat({er, pre_numtab.cols})
      i = i + 1
    end
  elseif er then
    er = tonumbers(er)
    vlst[#vlst+1] = nconcat({er, pre_numtab.cols})
    i = i + 1
  elseif nr then
    nr = tonumbers(nr)
    vlst[#vlst+1] = nconcat({pre_spectab.cols, nr})
    j = j + 1
  else
    break
  end
end
]]

local labels = {"Peaksuche", "Numerov", "exakt"}
local out = {}

local function write(...) out[#out+1] = string.format(...) end
--dump[1](vlst)

write("\\begin{tabular}{|")
for k = 1, #vlst[1] do
  write("c|")
end
write("}\\hline\n")
write("  \\rowcolor{lightgray}")
for k = 1, #vlst[1] do
  write(labels[k] or "")
  if k < #vlst[1] then
    write(" & ")
  end
end
write("\\\\ \\hline\n")
for k, row in ipairs(vlst) do
write("  ")
  for j, value in ipairs(row) do
    if value ~= "" then
      write(string.format("%11.10g", value))
    end
    if j < #row then
      write(" & ")
    end
  end
  write("\\\\ \\hline \n")
end

write("\\end{tabular}")

print(table.concat(out))

---