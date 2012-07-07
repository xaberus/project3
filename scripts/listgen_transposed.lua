
local dump = require "dump"

local tol = .01

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

assert(#args >= 2, "usage: tabgen numen.dat en1.dat en2.dat ...")

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

local samp = {}
local indices = {}

for k, v in ipairs(args) do
  samp[#samp+1] = read_table(v)
  indices[#indices+1] = 1
end

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
      local d = (min - ss[1]) / ss[1]
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

local out = {}

local function write(...) out[#out+1] = string.format(...) end
--dump[1](vlst)

local trans = {}
for k, row in ipairs(vlst) do
  local p = true
  for j, value in ipairs(row) do
    if value == "" then
      p = false
      break
    end
  end
  if p then
    trans[#trans+1] = row
  end
end

for k, row in ipairs(trans) do
  for j, value in ipairs(row) do
    if value ~= "" then
      write("%11.10g ", value)
    end
  end
  write("\n")
end
print(table.concat(out))

---
