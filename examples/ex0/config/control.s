# Example 0: Minimum settings
# ============================
# main parameters
# ============================
lnmlj   true
cutshortrange 10.0 
lverletL  true
trunc no
skindiff 1.5
# ============================
# md parameters 
# ============================
integrator nve-vv
temp 10.0
npas 0
dt .001
tauTberendsen .001
nprint 1
fprint 0
nequil 0
# ============================
# parameters for the force-field 
# ============================
massit 39.95
epslj   .010323576
sigmalj 3.405
