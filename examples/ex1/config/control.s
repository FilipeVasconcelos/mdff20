# Example 1: LJ fcc structure at low temperature
# ======================================
# main parameters
# ======================================
lnmlj          true    // switch on lj potential
cutshortrange  8.5    // lj short range cut off
lverletL       false 
# ======================================
# parameters for the molecular dynamics 
# ======================================
integrator     nve-vv  // NVE ensemble velocity verlet
npas           5000    // number of md steps
dt             .001   // time step    [ps]
temp           90.0   // temperature  [K]
nprint         5000   // printing info each nprint
fprint            1   // printing info to OSZIFF
cprint         5000   // printing config to CONTFF
nequil         1000    // equilibration for nequil steps
nequilT        1       // period of rescaling 
# ======================================
# parameters for the force-field 
# ======================================
epslj          .010323576 // lennard-jones potential in [eV] = 119.8 K
#sigmalj        3.82198327449341500178 
sigmalj 3.405      // lennard-jones potential in [A]
massit         39.95      // in atomic mass unit
