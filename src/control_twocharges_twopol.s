# ------------
#   global 
# ------------
Fposff rvf 
lverletL false 
cutshortrange 8.5 
cutlongrange  5.0 
skindiff  1.0 
lstatic true
lreduced true 
# ------------
#   field
# ------------
lcoulombic true
lautoES false 
alphaES 0.8
kES 11 11 11  
epsw 1e-10

lnmlj false 

# ------------
#    PIM 
# ------------
algo_pim scf
conv_tol_ind  1e-6    
algo_extrapolate_dipole aspc
extrapolate_order 4   
min_scf_pim_iter 3 
max_scf_pim_iter 3


#TYPE 1 
massit 1.0 
qit    1.0 
polit 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 

#TYPE 2 
massit 1.0 
qit   -1.0 
polit 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 

#TYPE 3 
massit 1.0
qit    0.0
polit 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 

#TYPE 4 
massit 1.0
qit    0.0
polit 0.1 0.0 0.0 0.0 0.1 0.0 0.0 0.0 0.1 
