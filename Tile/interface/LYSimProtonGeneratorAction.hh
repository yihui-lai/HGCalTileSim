#ifndef LYSimProtonGeneratorAction_h
#define LYSimProtonGeneratorAction_h

#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimDetectorConstruction.hh"
#else
#include "LYSimDetectorConstruction.hh"
#endif

#include "globals.hh"
#include "G4GeneralParticleSource.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UImessenger.hh"

class G4Event;
class LYSimProtonGeneratorMessenger;

class LYSimProtonGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
   LYSimProtonGeneratorAction();
   ~LYSimProtonGeneratorAction();

   void GeneratePrimaries( G4Event* );

   inline double
   GetBeamX() const { return _beamx; }
   inline void
   SetBeamX( const double x ){ _beamx = x; }
   inline double
   GetBeamY() const { return _beamy; }
   inline void
   SetBeamY( const double x ){ _beamy = x; }
   inline double
   GetWidth() const { return _width; }
   inline void
   SetWidth( const double x ){ _width = x; }

 private:
   G4GeneralParticleSource* particleSource;
   LYSimProtonGeneratorMessenger* messenger;

   double _beamx;
   double _beamy;
   double _width;

 };

 class G4UIdirectory;
 class G4UIcmdWithADoubleAndUnit;

 class LYSimProtonGeneratorMessenger : public G4UImessenger
 {
 public:
   LYSimProtonGeneratorMessenger( LYSimProtonGeneratorAction* );
   ~LYSimProtonGeneratorMessenger();

   void SetNewValue( G4UIcommand*, G4String );

 private:
   LYSimProtonGeneratorAction* genaction;

   G4UIdirectory* detDir;
   G4UIcmdWithADoubleAndUnit* SetBeamXCmd;
   G4UIcmdWithADoubleAndUnit* SetBeamYCmd;
   G4UIcmdWithADoubleAndUnit* SetWidthCmd;
 };


 #endif/*LYSimProtonGeneratorAction_h*/

