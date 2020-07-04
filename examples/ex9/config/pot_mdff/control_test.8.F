# ------------
#   Global 
# ------------
Fposff rvf
#lverletL false
lverletL true 
cutshortrange  7.2 
cutlongrange   7.2
skindiff  1.0 
lstatic false 
lreduced false 
lpstress true
# ------------
#   Field
# ------------
lcoulombic .true. 
lautoES true 
epsw     1e-12 
lbhmft  __BHM__
lbhmftd __BHM_DAMP__
massit 15.99994 # O
massit 28.084   # Si
massit 10.811   # B
qit    -2.0     # O     
qit     4.0     # Si 
qit     3.0     # B 
#--------------
#  BHMFTD
#--------------
lbhmftd true 
#         O-O             O-Si                O-B 
#         (eV)            (eV)                (eV)
Abhmftd 1.26277351E+04    1.31294927E+03      6.51712643E+02
Bbhmftd 5.0222184671669   3.1860795159        3.251153569410612
Cbhmftd   15.177181       1.1950537           0.0
Dbhmftd   82.256749       4.1831138           0.0
BDbhmftd  1.889           4.15739913          0.0
!INCLUDE <BMHFT.POT>  
# ------------
#    PIM 
# ------------
algo_pim scf
conv_tol_ind  1e-10    
algo_extrapolate_dipole aspc
extrapolate_order 22   
min_scf_pim_iter 10  
max_scf_pim_iter 100
!INCLUDE <PIM.POT>
# ----------------------------
# Polarizabilities
# ----------------------------
# O  type 
polit 0.93971813676926158987 0.0 0.0    0.0 0.93971813676926158987 0.0   0.0 0.0 0.93971813676926158987
# Si type
polit 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
# B type
polit 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 

# ----------------------------
#  Damping 
# ----------------------------
# Oxygen  O-O
lpoldamping 0 0 0 true 
pol_damp_b 0 0 0 5.58406251008641721012 
pol_damp_c 0 0 0 2.8904937083E+00
pol_damp_k 0 0 0 4

# Oxygen  O-Si
lpoldamping 0 0 1 true
pol_damp_b 0 0 1 3.66480
pol_damp_c 0 0 1 1.44589
pol_damp_k 0 0 1 4

# Oxygen  O-B
lpoldamping 0 0 2 true
pol_damp_b 0 0 2 4.48243215408077070621
pol_damp_c 0 0 2 1.3650000000E+00
pol_damp_k 0 0 2 4
