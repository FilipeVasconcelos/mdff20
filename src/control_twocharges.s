# ------------
#   global 
# ------------
Fposff 2
lverletL false 
cutshortrange 8.5 
cutlongrange  5.0 
skindiff  1.0 
lstatic true
lreduced true 
# ------------
#   field
# ------------
lcoul true
lautoES false 
qit 1.0000 -1.00000
alphaES 0.8
kES 11 11 11  
epsw 1e-10

lnmlj false 
massit 39.95

trunc no 
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
