include("sedt.jl")

using NetCDF

segmented = ncread(ARGS[1], "segmented")

tomo = compute(segmented)

(x, y, z) = size(tomo)
filename  = "tomo_float_sedt.nc"
varname   = "tomo_float"

nccreate(filename, varname,
         "$(varname)_xdim", x,
         "$(varname)_ydim", y,
         "$(varname)_zdim", z,
         t = NC_FLOAT)

ncwrite(tomo, filename, varname)
