# ------------
#   Global 
# ------------
Fposff rvf
lverletL true 
cutshortrange  5.0 
cutlongrange   5.0 
skindiff  1.0 
lstatic false 
lreduced false 
lpstress false 
# ------------
#   Field
# ------------
lcoulombic true 
lautoES true 
epsw 1e-6
massit 15.99994 # O
massit 28.084   # Si
qit    -2.0     # O     
qit     4.0     # Si 
# --------------
#   NMLJ
# --------------
#lnmlj true
#trunc no 
#massit 39.95
#epslj .010323576
#sigmalj 3.405
# --------------
#   BHMFTD
# --------------
lbhmftd true 
#         O-O        O-Si  
#         (eV)       (eV)  
Abhmftd 7902.1862   1311.5888
Bbhmftd    4.54668     3.18568
Cbhmftd   13.145548    0.0
Dbhmftd   71.279951    0.0
BDbhmftd   2.64562     0.0
# ------------
#   PIM 
# ------------
algo_pim scfkO
#algo_pim scf
omegakO 0.65
conv_tol_ind  1e-5    
algo_extrapolate_dipole aspc
extrapolate_order 4   
min_scf_pim_iter 4  
max_scf_pim_iter 100
# --------------------
#  Polarizability 
# --------------------
# O type 
polit 1.59150 0.0 0.0 0.0 1.59150 0.0 0.0 0.0 1.59150
#polit 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
# Si type
polit 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
# ------------
#  Damping O
# ------------
# O-O
lpoldamping 0 0 0 true 
pol_damp_b 0 0 0 4.74888
pol_damp_c 0 0 0 2.227
pol_damp_k 0 0 0 4
# O-Si
lpoldamping 0 0 1 true
pol_damp_b 0 0 1 3.66480
pol_damp_c 0 0 1 1.44589
pol_damp_k 0 0 1 4
# ------------
#    md  
# ------------
integrator nvt-nhcn
#integrator nve-vv
nhc_n          3
nhc_yosh_order 3
nhc_mults      2
timesca_thermo 1.0
temp 2500.0
npas 1000
nprint 100
fprint 10
cprint 10
dt .0005
tauTberendsen .0005
