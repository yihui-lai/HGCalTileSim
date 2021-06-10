#ifndef MATERIALS_HH
#define MATERIAL_HH

#include "G4Material.hh"

// list of functions for constructing non-standard material in Geant4
extern G4Material* Make_Bialkali();
extern G4Material* Make_Epoxy();
extern G4Material* Make_BC_630_grease();
extern G4Material* Make_EJ200();
extern G4Material* Make_SCSN81();
extern G4Material* Make_Resin();
extern void        Update_EJ200_AbsLength( G4Material*, const double x );
extern void        Update_flat_AbsLength( G4Material*, const double x );

extern G4Material* Make_Custom_Air();


#endif
