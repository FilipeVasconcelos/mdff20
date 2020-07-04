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
!INCLUDE <PIMD.POT>
