# ------------
#   Global 
# ------------
Fposff rvf
#lverletL false
lverletL true 
cutshortrange  __RCUT__ 
cutlongrange   __RCUT__
skindiff  1.0 
lstatic true 
lreduced false 
lpstress true
# ------------
#   Field
# ------------
lcoulombic __COUL__ 
lautoES true 
epsw     __EPSW__ 
lbhmft  __BHM__
lbhmftd __BHMD__
!INCLUDE <IONS.POT>
!INCLUDE <BHMFTD.POT>  
!INCLUDE <BHMFT.POT>  
# ------------
#    PIM 
# ------------
algo_pim scf
conv_tol_ind  1e-10    
algo_extrapolate_dipole aspc
extrapolate_order 2   
min_scf_pim_iter 10  
max_scf_pim_iter 100
!INCLUDE <PIM.POT>
!INCLUDE <PIMD.POT>
