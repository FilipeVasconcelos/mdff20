# ------------
#   global 
# ------------
Fposff 0
lverletL true 
cutshortrange 8.5 
skindiff  .5 
lstatic false
# ------------
#   field
# ------------
lnmlj true
trunc no 
massit 39.95
epslj .010323576
sigmalj 3.405
# ------------
#    md  
# ------------
# integrator nvt-nhcn
integrator nve-vv
nhc_yosh_order 3      
nhc_mults      3      
nhc_n          2 
timesca_thermo 10.0
temp 119.8 
npas 0
nprint 1000
fprint 10
nequil 1 
dt .004
tauTberendsen .004
