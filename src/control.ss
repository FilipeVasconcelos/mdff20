# Example 1: LJ fcc structure at low temperature
# ======================================
# main parameters
# ======================================
lnmlj          true    // switch on lj potential
cutshortrange  8.5F    // lj short range cut off
lverletL       true 
# ======================================
# parameters for the molecular dynamics 
# ======================================
integrator     nve-vv  // NVE ensemble velocity verlet
npas           100    // number of md steps
dt             .001F   // time step    [ps]
temp           90.0F   // temperature  [K]
nprint         10     // printing info each nprint
fprint         1       // printing info to OSZIFF
nequil         1000    // equilibration for nequil steps
nequil_period  1       // period of rescaling 
# ======================================
# parameters for the force-field 
# ======================================
epslj     .010323576F // lennard-jones potential in [eV] = 119.8 K
sigmalj   3.405F      // lennard-jones potential in [A]
massit    39.95F      // in atomic mass unit
