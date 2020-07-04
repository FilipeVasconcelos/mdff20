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
!INCLUDE <BMHFTD.POT>  
#--------------
#  BHMFTD
#--------------
#         O-O             O-Si                O-B 
#         (eV)            (eV)                (eV)
Abhmftd 1.26277351E+04    1.31294927E+03      6.51712643E+02
Bbhmftd 5.0222184671669   3.1860795159        3.251153569410612
Cbhmftd   15.177181       1.1950537           0.0
Dbhmftd   82.256749       4.1831138           0.0
BDbhmftd   0.0            0.0                 0.0
# ------------
#    PIM 
# ------------
algo_pim scf
conv_tol_ind  1e-10    
algo_extrapolate_dipole aspc
extrapolate_order 22   
min_scf_pim_iter 10  
max_scf_pim_iter 100
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
lpoldamping 0 0 0 false 
# Oxygen  O-Si
lpoldamping 0 0 1 false 
# Oxygen  O-B
lpoldamping 0 0 2  false
!INCLUDE <PIMD.POT>
