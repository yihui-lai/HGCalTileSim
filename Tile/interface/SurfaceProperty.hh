#ifndef SURFACE_PROPERTY_HH
#define SURFACE_PROPERTY_HH

#include "G4OpticalSurface.hh"

extern G4OpticalSurface* MakeS_TyvekCrystal();
extern G4OpticalSurface* MakeS_ESR();
extern G4OpticalSurface* MakeS_IdealTyvekCrystal();
extern G4OpticalSurface* MakeS_UnifiedTyvek();
extern G4OpticalSurface* MakeS_IdealTyvek();
extern G4OpticalSurface* MakeS_Polished();
extern G4OpticalSurface* MakeS_IdealPolished();
extern G4OpticalSurface* MakeS_Mirror();
extern G4OpticalSurface* MakeS_IdealMirror();
extern G4OpticalSurface* MakeS_IdealWhiteSurface();
extern G4OpticalSurface* MakeS_Absorbing();
extern G4OpticalSurface* MakeS_RoughMirror();
extern G4OpticalSurface* MakeS_RoughInterface( double );
extern G4OpticalSurface* MakeS_SiPM();
extern G4OpticalSurface* MakeS_PCBSurface();
extern G4OpticalSurface* MakeS_AirAir();

//wls

extern G4OpticalSurface* MakeS_TiO2Surface();


#endif