# ------------
#   global 
# ------------
Fposff rvf
#lverletL false
lverletL true 
cutshortrange  5.5 
cutlongrange   5.0 
skindiff  1.0 
lstatic false 
lreduced false 
lpstress false 
# ------------
#   field
# ------------
lcoulombic true 
lautoES true 
epsw 1e-6

#--------------
# NMLJ
#--------------
#lnmlj true
trunc no 
#massit 39.95
#epslj .010323576
#sigmalj 3.405
#--------------
#  BHMFTD
#--------------
lbhmft true 
massit 15.99994 # O
massit 28.084   # Si
qit -2.0 #O     
qit  4.0 #Si 

#         O-O        O-Si  
#         (eV)       (eV)  
Abhmftd 7902.1862   1311.5888
Bbhmftd    4.54668     3.18568
Cbhmftd   13.145548    0.0
Dbhmftd   71.279951    0.0
# ------------
#    PIM 
# ------------
algo_pim scf
conv_tol_ind  1e-6    
algo_extrapolate_dipole aspc
extrapolate_order 4   

# Oxygen type 
lpolar true
polit 1.59150 0.0 0.0 0.0 1.59150 0.0 0.0 0.0 1.59150
lpolar false 

# ------------
#    md  
# ------------
integrator nvt-nhcn
#integrator nve-vv
nhc_n          4 
nhc_yosh_order 5      
nhc_mults      4      
timesca_thermo 1.0 
temp 2500.0 
npas 0 
nprint 1000
fprint 10
cprint 100000
dt .0005
tauTberendsen .0005
