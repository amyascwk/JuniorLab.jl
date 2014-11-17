using Base.Test
using JLab
import PyPlot.plt

#Sample data
y = [0,0,4,4,4,4,6,8,10,11,10,8,6,4,2,0,0,0,-2,-4,-6,-4,-6,-4,-2,0,2,4,6,8,10,4,6,4,2,0,0,-2,-4,-4]
peakloc = [10,22,31,33]
lbnds_h = [6,21,23,32]
lbnds_m = [2,21,23,32]
rbnds_h = [16,23,32,36]
rbnds_m = [21,23,32,39]

#  plt.plot(y)

# #############################################################################
#Testing naive peak finding
println("########## Testing naive peak finding ##########")

#Testing hills
println("Testing hills ...")
pks_h = gethills(y)
@test pks_h.summits == peakloc
@test pks_h.lbnds == lbnds_h
@test pks_h.rbnds == rbnds_h

#Testing mountains
println("Testing mountains ...")
pks_m = getmountains(y)
@test pks_m.summits == peakloc
@test pks_m.lbnds == lbnds_m
@test pks_m.rbnds == rbnds_m

#Testing threshhold
println("Testing threshhold ...")
pks_h = gethills(y,thresh=2)
@test pks_h.summits == [10]
@test pks_h.lbnds == [2]
@test pks_h.rbnds == [32]
pks_h = gethills(y,thresh=4)
@test pks_h.summits == []

pks_m = getmountains(y,thresh=2)
@test pks_m.summits == [10]
@test pks_m.lbnds == [2]
@test pks_m.rbnds == [32]
pks_m = getmountains(y,thresh=4)
@test pks_m.summits == []

println("Functions \"getmountains\" and \"gethills\" passed tests.")

