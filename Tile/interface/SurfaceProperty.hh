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
extern G4OpticalSurface* MakeS_DarkSurface();
extern G4OpticalSurface* MakeS_Absorbing();
extern G4OpticalSurface* MakeS_Rough();
extern G4OpticalSurface* MakeS_SiPM_S13360();
extern G4OpticalSurface* MakeS_SiPM_S14160();
extern G4OpticalSurface* MakeS_PCBSurface();
extern G4OpticalSurface* MakeS_AirAir();
extern G4OpticalSurface* MakeS_dimple();
extern void  Update_Dimple_SigmaAlpha( G4OpticalSurface*, const double x );
#endif
