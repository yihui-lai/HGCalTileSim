#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimDetectorConstruction.hh"
#include "HGCalTileSim/Tile/interface/LYSimDetectorMessenger.hh"
#include "HGCalTileSim/Tile/interface/Materials.hh"
#include "HGCalTileSim/Tile/interface/SurfaceProperty.hh"
#else
#include "LYSimDetectorConstruction.hh"
#include "LYSimDetectorMessenger.hh"
#include "Materials.hh"
#include "SurfaceProperty.hh"
#endif

#include <math.h>
#include <vector>

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4GeometryManager.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4RegionStore.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4MultiUnion.hh"
#include "G4UnionSolid.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"
#include "TF1.h"
#include "G4Polycone.hh"

using std::cos;
using std::sin;
using std::tan;
using std::atan;
using std::exp;
using namespace CLHEP;

LYSimPMTSD* LYSimDetectorConstruction::fPMTSD = NULL;

LYSimDetectorConstruction::LYSimDetectorConstruction()
  : G4VUserDetectorConstruction()
{
  fdetectorMessenger = new LYSimDetectorMessenger( this );

  _tilex        = 30*mm;//29.65*mm;
  _tiley        = 30*mm;//29.65*mm;
  _tilez        = 3.0*mm;//3.0*mm;
  _tile_x1      = 0.0*mm;
  _tile_x2      = 0.0*mm;

  wrapgap       = 0.3*mm;
  wrapthickness = 0.1*mm;

  tilegap = 0.1*mm;
  tile_number = 1;
  is_ESR = true;
  material = 1;
  _dimplesa = 1.3;

  _absmult      = 1;
  _wrap_reflect = 0.985;

//this is the 1.3mm SiPM
  _sipm_deadwidth  = 0.2*mm;
  _sipm_x          = 1.3*mm;
  _sipm_y          = 1.3*mm;
  _sipm_z          = 0.35*mm;
  _sipm_rimwidth   = 0.4*mm;
  _sipm_glasswidth = 0.3*mm;
  _sipm_standz     = 0.2*mm;
  _gap_pcb_wrap    = 0*mm;

/*
//this is the 1.4mm SiPM
  _sipm_deadwidth  = 0.2*mm;
  _sipm_x          = 1.4*mm;
  _sipm_y          = 1.4*mm;
  _sipm_z          = 0.4*mm;
  _sipm_rimwidth   = 0.3*mm;
  _sipm_glasswidth = 0.1*mm;
  _sipm_standz     = 0.3*mm;

//this is the 2.0mm SiPM
  _sipm_x          = 2.0*mm;
  _sipm_y          = 2.0*mm;
  _sipm_rimwidth   = 0.1*mm;
*/
/*
//this is the 3.0mm SiPM
  _sipm_x          = 3.0*mm;
  _sipm_y          = 3.0*mm;
  _sipm_rimwidth   = 0.1*mm;
*/

  // Default Dimple settings
  _dimple_type   = SPHERICAL;// 0: Normal, 1: Pyramid, 2: Parabolic
  _dimple_indent = 1.6*mm;
  _dimple_radius = 2.65*mm;

  // Default Hole settings
  _pcb_radius       = 3.1;
  _pcb_reflectivity = 0.5;

  // Defining material list.
  fBialkali = Make_Bialkali();
  fEpoxy    = Make_Epoxy();
  fAir      = Make_Custom_Air();
  if(material==1){ fTileMaterial = Make_EJ200();}else{fTileMaterial = Make_SCSN81();}
 
  fResin    = Make_Resin();
  fBC_630_grease = Make_BC_630_grease();
  //SetTileAbsMult( _absmult );
  // Defining surface list.
  fTyvekSurface           = MakeS_TyvekCrystal();
  fESROpSurface           = MakeS_Rough();
  fPolishedOpSurface      = MakeS_Polished();
  fIdealPolishedOpSurface = MakeS_IdealPolished();
  fIdealWhiteOpSurface    = MakeS_IdealWhiteSurface();
  fDarkOpSurface    = MakeS_DarkSurface();
  fSiPMSurface            = MakeS_SiPM_S13360();
  fPCBSurface             = MakeS_PCBSurface();
  fGrease                 = MakeS_Absorbing();
  fdimpleSurface          = MakeS_IdealPolished();
  //SetWrapReflect( _wrap_reflect );

}

void
LYSimDetectorConstruction::UpdateGeometry()
{
  
  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();

  G4RunManager::GetRunManager()->DefineWorldVolume( Construct() );
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

LYSimDetectorConstruction::~LYSimDetectorConstruction()
{
  if( fdetectorMessenger ){ delete fdetectorMessenger; }
}

G4VPhysicalVolume*
LYSimDetectorConstruction::Construct()
{
  static const bool checkOverlaps = true;

std::cout<<"############## tile: ("<<_tilex<<" "<<_tiley<<" "<<_tilez<<")"<<std::endl;
std::cout<<"############## wrap gap: "<<wrapgap<<" thickness: "<<wrapthickness<<std::endl;
std::cout<<"############## tilegap: "<<tilegap<<" number: "<<tile_number<<std::endl;
std::cout<<"############## material: "<<material<<" is esr:"<<is_ESR<<" _dimplesa: "<<_dimplesa<<std::endl;
std::cout<<"############## absmul: "<<_absmult<<" wrap_ref:"<<_wrap_reflect<<std::endl;
std::cout<<"############## pde: "<<sipm_pde<<std::endl;
if(material==0){fTileMaterial = Make_SCSN81();}
if(sipm_pde==1){fSiPMSurface            = MakeS_SiPM_S13360();}else{fSiPMSurface            = MakeS_SiPM_S14160();}
SetTileAbsMult( _absmult ); std::cout<<"set TileAbsMult:"<<_absmult<<std::endl;
SetWrapReflect( _wrap_reflect ); std::cout<<"set wrap reflect:"<<_wrap_reflect<<std::endl;

  ///////////////////////////////////////////////////////////////////////////////
  // World volume
  ///////////////////////////////////////////////////////////////////////////////
  G4Box* solidWorld = new G4Box( "World", WorldHalfX()
                               , WorldHalfY(), WorldHalfZ() );

  G4LogicalVolume* logicWorld = new G4LogicalVolume( solidWorld, fAir, "World" );

  G4VPhysicalVolume* physWorld = new G4PVPlacement( 0
                                                  , G4ThreeVector( 0, 0, 0 )
                                                  , logicWorld
                                                  , "World"
                                                  , 0
                                                  , false
                                                  , 0
                                                  , checkOverlaps );
  ///////////////////////////////////////////////////////////////////////////////
  // wrapping
  ///////////////////////////////////////////////////////////////////////////////
  //incase the wrapper is too small
  G4VSolid* solidWrap_pre = ConstructHollowWrapSolid();

  G4ThreeVector SiPMOffset_pre( -0.5*0.53*mm, 0,  - (0.5*_sipm_z + _sipm_standz + wrapthickness) + 0.5*_tilez + wrapgap + 2*wrapthickness +_gap_pcb_wrap + wrapgap - 0.1*mm);

  G4Box* solidSiPMResinOuter_pre
    = new G4Box( "solidSiPMResinOuter_pre"
               , 0.5*_sipm_x + _sipm_rimwidth + 0.5*0.53*mm
               , 0.5*_sipm_y + _sipm_rimwidth
               , 3*_sipm_z + 3*_sipm_glasswidth );
  G4VSolid* solidWrap = new G4SubtractionSolid( "WrapBox"
                                            , solidWrap_pre, solidSiPMResinOuter_pre , NULL, SiPMOffset_pre );


  // The matrial of the wrap isn't as important as the surface
  G4LogicalVolume* logicWrap = new G4LogicalVolume( solidWrap, fEpoxy,  "Wrap" );

  const G4ThreeVector tile_offset( 0, 0
                                , -0.5*(_tilez*(tile_number-1) + tilegap*(tile_number-1)) + wrapgap - 0.1*mm);

  G4VPhysicalVolume* physWrap = new G4PVPlacement( 0
                                                 , tile_offset
                                                 , logicWrap
                                                 , "Wrap"
                                                 , logicWorld
                                                 , false
                                                 , 0
                                                 , checkOverlaps );
  // Exposed PCB Back plane
  const G4ThreeVector pcb_offset( 0, 0
                                , 0.5*_tilez + wrapgap + 2*wrapthickness +_gap_pcb_wrap + wrapgap - 0.1*mm );

  G4VSolid* solidPCB = new G4Tubs( "PCBSolid"
                                 , 0, _pcb_radius + 1*mm
                                 , wrapthickness
                                 , 0, 2*pi );
  G4LogicalVolume* logicPCB = new G4LogicalVolume( solidPCB, fEpoxy, "PCB" );
/*
  G4VPhysicalVolume* physPCB = new G4PVPlacement( 0, pcb_offset
                                                , logicPCB, "PCB"
                                                , logicWorld, false
                                                , 0, checkOverlaps  );
*/
  ///////////////////////////////////////////////////////////////////////////////
  //  Add two proton position trackers
  ///////////////////////////////////////////////////////////////////////////////
  G4Box* solidTracker = new G4Box( "Tracker", 1.5*_tilex, 1.5*_tiley, 0.1*mm  );
  G4LogicalVolume* logicTracker = new G4LogicalVolume( solidTracker, fEpoxy, "Tracker" );
  G4VPhysicalVolume* physTrackerFront = new G4PVPlacement( 0
                                                  , pcb_offset + G4ThreeVector( 0, 0, 10*mm )
                                                  , logicTracker
                                                  , "TrackerFront"
                                                  , logicWorld
                                                  , false
                                                  , 0
                                                  , checkOverlaps );
  G4VPhysicalVolume* physTrackerBack = new G4PVPlacement( 0
                                                  , G4ThreeVector( 0, 0, -(_tilez*(tile_number-1) + tilegap*(tile_number-1))-10*mm )
                                                  , logicTracker
                                                  , "TrackerBack"
                                                  , logicWorld
                                                  , false
                                                  , 0
                                                  , checkOverlaps );


  ///////////////////////////////////////////////////////////////////////////////
  // Subtracted Dimple Version (dimple sub from tile, WWW = mothervolume of both)
  ///////////////////////////////////////////////////////////////////////////////
  G4VSolid* solidTile = ConstructTrapazoidSolid( "TileTrap"
                                               , _tilex
                                               , _tiley
                                               , _tilez
                                               , _tile_x1
                                               , _tile_x2 );

  G4VSolid* solidTileSubtract = new G4Sphere( "TileSub", 0, GetDimpleSizeRadius()
                                            , 0, 2* pi, 0, pi );
  G4VSolid* hollowedTile = new G4SubtractionSolid( "TileSolid"
                                                 , solidTile, solidTileSubtract
                                                 , 0, CalcDimpleOffset() );

  //replace the above code with the following to construct any shape
  /*
  G4VSolid* solidTileSubtract = ConstructArbiDimpleSolid();
  G4VSolid* hollowedTile = new G4SubtractionSolid( "TileSolid"
                                                 , solidTile, solidTileSubtract
                                                 , 0, G4ThreeVector( 0, 0, 0.5*_tilez ) );
  */


  G4LogicalVolume* logicTile = new G4LogicalVolume( hollowedTile
                                                  , fTileMaterial, "TileLogic" );

  G4VPhysicalVolume* physTile = new G4PVPlacement( 0
                                                 , G4ThreeVector( 0, 0, 0 )
                                                 , logicTile
                                                 , "TilePhysic"
                                                 , logicWorld
                                                 , false
                                                 , 0
                                                 , checkOverlaps );

  ///////////////////////////////////////////////////////////////////////////////
  // put more tile layers and grease
  ///////////////////////////////////////////////////////////////////////////////

  G4LogicalVolume* mlogicTile= new G4LogicalVolume( solidTile
                                               , fTileMaterial, "mTileLogic" );

  std::vector<G4VPhysicalVolume*> mphysTileList;

  G4VSolid* msolidGrease = ConstructTrapazoidSolid( "GreaseTrap"
                                               , _tilex
                                               , _tiley
                                               , tilegap
                                               , _tile_x1
                                               , _tile_x2 );
  G4LogicalVolume* mlogicGrease= new G4LogicalVolume( msolidGrease
                                               , fBC_630_grease, "GreaseLogic" );

  G4LogicalSkinSurface* mGreaseAbsorb = new G4LogicalSkinSurface(
                               "mAbsorbGreaseSurface"
                              , mlogicGrease
                              , fGrease);

  std::vector<G4VPhysicalVolume*> mphysGreaseList;

  for(int i=1; i<tile_number; i++){
      mphysGreaseList.push_back( 
        new G4PVPlacement(0, G4ThreeVector( 0, 0, -1*(_tilez+tilegap)*i+0.5*(_tilez+tilegap) )
                                               , mlogicGrease
                                               , "GreasePhysic"
                                               , logicWorld
                                               , false
                                               , 0
                                               , checkOverlaps)); 

      mphysTileList.push_back( new G4PVPlacement(0
                                               , G4ThreeVector( 0, 0, -1*(_tilez+tilegap)*i )
                                               , mlogicTile
                                               , "mTilePhysic"
                                               , logicWorld
                                               , false
                                               , 0
                                               , checkOverlaps ));
  }



  ///////////////////////////////////////////////////////////////////////////////
  // this is the start of 1.3 mm SiPM
  ///////////////////////////////////////////////////////////////////////////////

  double extra_l = 0.53*mm;
  double thin_surrounding = 0.01*mm;

  G4Box* solidSiPMactive = new G4Box( "SiPMDead"
                                  , 0.5*_sipm_x, 0.5*_sipm_y
                                  , _sipm_z );

  G4Box* solidSiPMInnerBox = new G4Box( "SiPMInnerBox"
                                      , 0.5*_sipm_x+thin_surrounding, 0.5*_sipm_y+thin_surrounding
                                      ,  0.5*_sipm_z);


  G4VSolid* solidSiPMCase
    = new G4SubtractionSolid( "SiPMCase"
                            , solidSiPMInnerBox, solidSiPMactive
                            , 0
                            , G4ThreeVector( 0, 0, 0 ) );

  G4VSolid* solidSiPMInner
    = new G4IntersectionSolid( "SiPMInner"
                             , solidSiPMInnerBox, solidSiPMactive
                             , 0
                             , G4ThreeVector( 0, 0, 0 ) );

  G4Box* solidSiPMResinOuter
    = new G4Box( "SiPMResinOuter"
               , 0.5*_sipm_x + _sipm_rimwidth + 0.5*extra_l
               , 0.5*_sipm_y + _sipm_rimwidth
               , 0.5*_sipm_z + 0.5*_sipm_glasswidth );

  G4VSolid* solidSiPMResin
    = new G4SubtractionSolid( "SiPMResin"
                            , solidSiPMResinOuter, solidSiPMInnerBox
                            ,  0
                            , G4ThreeVector( 0.5*extra_l, 0, 0.5*_sipm_glasswidth ) );
  G4Box* solidSiPMStand
    = new G4Box( "SiPMStand"
               , 0.5*_sipm_x+_sipm_rimwidth + 0.5*extra_l
               , 0.5*_sipm_y+_sipm_rimwidth
               , 0.5*_sipm_standz );


  G4LogicalVolume* logicSiPM = new G4LogicalVolume( solidSiPMInner
                                                  , fBialkali,  "SiPM" );

  G4LogicalVolume* logicSiPMCase = new G4LogicalVolume( solidSiPMCase
                                                      , fEpoxy, "SiPMBack" );

  G4LogicalVolume* logicSiPMResin = new G4LogicalVolume( solidSiPMResin
                                                       , fResin, "SiPMResin" );

  G4LogicalVolume* logicSiPMStand = new G4LogicalVolume( solidSiPMStand
                                                       , fEpoxy, "SiPMStand" );

  G4ThreeVector SiPMOffset( 0, 0,  - (0.5*_sipm_z + _sipm_standz + wrapthickness) );
  SiPMOffset +=pcb_offset;
  G4ThreeVector ResinOffset( -0.5*extra_l, 0,  - (0.5*_sipm_glasswidth + 0.5*_sipm_z + _sipm_standz + wrapthickness) );
  ResinOffset+=pcb_offset;
  G4ThreeVector StandOffset( -0.5*extra_l, 0,  - (0.5*_sipm_standz + wrapthickness) );
  StandOffset+=pcb_offset;

  G4VPhysicalVolume* physSiPMStand = new G4PVPlacement( 0
                                                      , StandOffset
                                                      , logicSiPMStand
                                                      , "SiPMStand"
                                                      , logicWorld
                                                      , false
                                                      , 0
                                                      , checkOverlaps );

  G4VPhysicalVolume* physSiPMCase = new G4PVPlacement( 0
                                                     , SiPMOffset
                                                     , logicSiPMCase
                                                     , "Case"
                                                     , logicWorld
                                                     , false
                                                     , 0
                                                     , checkOverlaps );

  G4VPhysicalVolume* physSiPMResin = new G4PVPlacement( 0
                                                      , ResinOffset
                                                      , logicSiPMResin
                                                      , "SiPMResin"
                                                      , logicWorld
                                                      , false
                                                      , 0
                                                      , checkOverlaps  );

  G4VPhysicalVolume* physSiPM = new G4PVPlacement( 0
                                                 , SiPMOffset
                                                 , logicSiPM
                                                 , "SiPM"
                                                 , logicWorld
                                                 , false
                                                 , 0
                                                 , checkOverlaps );




///////////////////////////////////////////////////////////////////////////////
// this is the 1.4 mm SiPM
///////////////////////////////////////////////////////////////////////////////
/*
  G4Box* solidSiPMDead = new G4Box( "SiPMDead"
                                  , 0.5*_sipm_deadwidth, 0.5*_sipm_deadwidth
                                  , _sipm_z );

  G4Box* solidSiPMInnerBox = new G4Box( "SiPMInnerBox"
                                  , 0.5*_sipm_x, 0.5*_sipm_y,  0.8*_sipm_z );

  G4Box* solidSiPMOuter = new G4Box( "SiPMOuter"
                                  , 0.5*_sipm_x + _sipm_rimwidth
                                  , 0.5*_sipm_y + _sipm_rimwidth
                                  , 0.5*_sipm_z );
  G4Box* solidSiPMStand = new G4Box( "SiPMStand"
                                  , 0.5*_sipm_x+_sipm_rimwidth + _sipm_glasswidth
                                  , 0.5*_sipm_y+_sipm_rimwidth + _sipm_glasswidth
                                  , 0.5*_sipm_standz );

  G4Box* solidSiPMResinOuter
    = new G4Box( "SiPMResinOuter"
               , 0.5*_sipm_x + _sipm_rimwidth + _sipm_glasswidth
               , 0.5*_sipm_y + _sipm_rimwidth + _sipm_glasswidth
               , 0.5*_sipm_z + _sipm_glasswidth );
  G4VSolid* solidSiPMSubtract
    = new G4SubtractionSolid( "SiPMSubtract"
                            ,  solidSiPMInnerBox, solidSiPMDead
                            ,   0, G4ThreeVector( 0, 0, 0 ) );
  G4VSolid* solidSiPMCase
    = new G4SubtractionSolid( "SiPMCase"
                            , solidSiPMOuter, solidSiPMSubtract
                            , 0
                            , G4ThreeVector( 0, 0, -0.65 * _sipm_z ) );

  G4VSolid* solidSiPMInner
    = new G4IntersectionSolid( "SiPMInner"
                             , solidSiPMOuter, solidSiPMSubtract
                             , 0
                             , G4ThreeVector( 0, 0, -0.65*_sipm_z ) );

  G4VSolid* solidSiPMResin
    = new G4SubtractionSolid( "SiPMResin"
                            , solidSiPMResinOuter, solidSiPMOuter
                            ,  0
                            , G4ThreeVector( 0, 0, _sipm_glasswidth ) );
  G4LogicalVolume* logicSiPM = new G4LogicalVolume( solidSiPMInner
                                                  , fBialkali,  "SiPM" );

  G4LogicalVolume* logicSiPMCase = new G4LogicalVolume( solidSiPMCase
                                                      , fEpoxy, "SiPMBack" );

  G4LogicalVolume* logicSiPMResin = new G4LogicalVolume( solidSiPMResin
                                                       , fResin, "SiPMResin" );

  G4LogicalVolume* logicSiPMStand = new G4LogicalVolume( solidSiPMStand
                                                       , fEpoxy, "SiPMStand" );
  
  G4ThreeVector SiPMOffset( 0, 0,  - (0.5*_sipm_z + _sipm_standz + wrapthickness) );
  SiPMOffset +=pcb_offset;
  G4ThreeVector ResinOffset( 0, 0,  - (_sipm_glasswidth + 0.5*_sipm_z + _sipm_standz + wrapthickness) );
  ResinOffset+=pcb_offset;
  G4ThreeVector StandOffset( 0, 0,  - (0.5*_sipm_standz + wrapthickness) );
  StandOffset+=pcb_offset;
  
  
  G4VPhysicalVolume* physSiPMStand = new G4PVPlacement( 0
                                                      , StandOffset
                                                      , logicSiPMStand
                                                      , "SiPMStand"
                                                      , logicWorld
                                                      , false
                                                      , 0
                                                      , checkOverlaps );

  G4VPhysicalVolume* physSiPMCase = new G4PVPlacement( 0
                                                     , SiPMOffset
                                                     , logicSiPMCase
                                                     , "Case"
                                                     , logicWorld
                                                     , false
                                                     , 0
                                                     , checkOverlaps );

  G4VPhysicalVolume* physSiPMResin = new G4PVPlacement( 0
                                                      , ResinOffset
                                                      , logicSiPMResin
                                                      , "SiPMResin"
                                                      , logicWorld
                                                      , false
                                                      , 0
                                                      , checkOverlaps  );

G4VPhysicalVolume* physSiPM = new G4PVPlacement( 0
                                                 , SiPMOffset
                                                 , logicSiPM
                                                 , "SiPM"
                                                 , logicWorld
                                                 , false
                                                 , 0
                                                 , checkOverlaps );


*/


  ///////////////////////////////////////////////////////////////////////////////
  // Defining surfaces
  ///////////////////////////////////////////////////////////////////////////////
  if (is_ESR){
        G4LogicalBorderSurface*  WrapAirSurface = new G4LogicalBorderSurface( "WrapAirSurface"
                              , physWorld
                              , physWrap
                              , fESROpSurface );
  }else{
        G4LogicalBorderSurface*  WrapAirSurface = new G4LogicalBorderSurface( "WrapAirSurface"
                              , physWorld
                              , physWrap
                              , fTyvekSurface );
  }

  G4LogicalBorderSurface* TileSurface =
    new G4LogicalBorderSurface( "TileSurface"
                              , physTile
                              , physWorld
                              , fIdealPolishedOpSurface );
  std::vector<G4LogicalBorderSurface*> mTileSurfaceList;

  for(int i=1; i<tile_number; i++){
    mTileSurfaceList.push_back( new G4LogicalBorderSurface(
                               "mTileSurface"
                              , mphysTileList.at(i-1)
                              , physWorld
                              , fIdealPolishedOpSurface) );    
  }
  //use white surface for 1.4mm, dark surface for 1.3mm
  G4LogicalSkinSurface* CaseSurface
    = new G4LogicalSkinSurface( "SiPMCaseSurface"
                              , logicSiPMCase
                              , fDarkOpSurface );
//                              , fIdealWhiteOpSurface );
  G4LogicalSkinSurface* StandSurface
    = new G4LogicalSkinSurface( "SiPMStandSurface"
                              , logicSiPMStand
                              , fIdealWhiteOpSurface );

  G4LogicalSkinSurface* SiPMSurface
    = new G4LogicalSkinSurface( "SiPMSurface", logicSiPM, fSiPMSurface );
  G4LogicalSkinSurface* PCBSurface
    = new G4LogicalSkinSurface( "PCBSurface", logicPCB, fPCBSurface );

  // Setting the sensitive detector
  if( !fPMTSD ){
    fPMTSD = new LYSimPMTSD( "/LYSimPMT" );
    G4SDManager* sdman = G4SDManager::GetSDMpointer();
    sdman->AddNewDetector( fPMTSD );
  }
  logicSiPM->SetSensitiveDetector( fPMTSD );

  // Visual attributes
  logicWorld->SetVisAttributes( G4VisAttributes::Invisible );

  G4VisAttributes* SiPMVisAtt = new G4VisAttributes( G4Colour( 0, 0, 0 ) );
  SiPMVisAtt->SetForceSolid( true );
  SiPMVisAtt->SetVisibility( true );
  logicSiPM->SetVisAttributes( SiPMVisAtt );

  G4VisAttributes* CaseVisAtt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0.8 ) );
  CaseVisAtt->SetForceSolid( true );
  CaseVisAtt->SetVisibility( true );
  logicSiPMCase->SetVisAttributes( CaseVisAtt );
  logicSiPMStand->SetVisAttributes( CaseVisAtt );

  G4VisAttributes* ResinVisAtt = new G4VisAttributes( G4Colour( 0., 1., 1. ) );
  ResinVisAtt->SetForceWireframe( true );
  ResinVisAtt->SetVisibility( true );
  logicSiPMResin->SetVisAttributes( ResinVisAtt );

  G4VisAttributes* TileVisAtt = new G4VisAttributes( G4Colour( 1., 1., 0. ) );
  TileVisAtt->SetForceWireframe( true );
  TileVisAtt->SetVisibility( true );
  logicTile->SetVisAttributes( TileVisAtt );

  mlogicTile->SetVisAttributes( TileVisAtt );

  G4VisAttributes* GreaseVisAtt = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0 ) );
  GreaseVisAtt->SetVisibility( true );
  mlogicGrease->SetVisAttributes( GreaseVisAtt );

  G4VisAttributes* WrapVisAtt = new G4VisAttributes( G4Colour( 0.5, 1., 0.5 ) );
  WrapVisAtt->SetForceWireframe( true );
  WrapVisAtt->SetVisibility( true );
  logicWrap->SetVisAttributes( WrapVisAtt );

  G4VisAttributes* TrackerVisAtt = new G4VisAttributes( G4Colour(1., 1., 0.) );
  TrackerVisAtt->SetForceWireframe( true );
  TrackerVisAtt->SetVisibility( true );
  logicTracker->SetVisAttributes( TrackerVisAtt );
 
  G4VisAttributes* PCBVisAtt = new G4VisAttributes( G4Colour( 0.0, 0.4, 0.1 ) );
  PCBVisAtt->SetForceSolid( true );
  PCBVisAtt->SetVisibility( true );
  logicPCB->SetVisAttributes( PCBVisAtt );

  return physWorld;
}

G4VSolid*
LYSimDetectorConstruction::ConstructTrapazoidSolid(
  const G4String& name,
  double          x,
  double          y,
  double          z,
  double          indent_x1,
  double          indent_x2 ) const
{
  const G4ThreeVector corners[8] = {
    G4ThreeVector( -x/2,           -y/2, -z/2 ),
    G4ThreeVector( x/2,            -y/2, -z/2 ),
    G4ThreeVector( -x/2+indent_x2, y/2,  -z/2 ),
    G4ThreeVector( x/2-indent_x1,  y/2,  -z/2 ),
    G4ThreeVector( -x/2,           -y/2, z/2 ),
    G4ThreeVector( x/2,            -y/2, z/2 ),
    G4ThreeVector( -x/2+indent_x2, y/2,  z/2 ),
    G4ThreeVector( x/2-indent_x1,  y/2,  z/2 )
  };

  return new G4Trap( name, corners );
}

G4VSolid*
LYSimDetectorConstruction::ConstructHollowWrapSolid() const
{
  G4VSolid* wrapOuter
    = ConstructTrapazoidSolid( "WrapOuter"
                             , _tilex + 2*0.1*mm + 2*wrapthickness
                             , _tiley + 2*0.1*mm + 2*wrapthickness
                             , _tilez*tile_number + tilegap*(tile_number-1) + 2*wrapgap + 2*wrapthickness
                             , 0, 0 );
  G4VSolid* wrapInner
    = ConstructTrapazoidSolid( "WrapInner"
                             , _tilex + 2*0.1*mm
                             , _tiley + 2*0.1*mm
                             , _tilez*tile_number + tilegap*(tile_number-1) + 2*wrapgap
                             , 0, 0 );
  G4VSolid* wrapbox = new G4SubtractionSolid( "WrapBox"
                                            , wrapOuter, wrapInner );
  G4VSolid* wraphole = new G4Tubs( "WrapHole"
                                 , 0, _pcb_radius
                                 , _tilez
                                 , 0, 2*pi );

  const G4ThreeVector offset( 0, 0, 0.5*(_tilez*tile_number + tilegap*(tile_number-1)) + wrapgap + 0.5*wrapthickness );
  G4EllipticalTube* wrapoval = new G4EllipticalTube("WrapHole", 1.9*mm, 1.65*mm, _tilez);

  if(_pcb_radius <= 1.8*mm ){ return new G4SubtractionSolid( "WrapSolid", wrapbox, wrapoval, NULL, offset );}
  else{ return new G4SubtractionSolid( "WrapSolid", wrapbox, wraphole, NULL, offset );}

}


G4VSolid*
LYSimDetectorConstruction::ConstructUnEvenHollowWrapSolid() const
{
  G4VSolid* wrapOuter
    = ConstructTrapazoidSolid( "WrapOuter"
                             , _tilex + wrapgap + 2*wrapthickness
                             , _tiley + wrapgap + 2*wrapthickness
                             , _tilez*tile_number + tilegap*(tile_number-1) + wrapgap + 2*wrapthickness
                             , 0, 0 );
  G4VSolid* wrapInner
    = ConstructTrapazoidSolid( "WrapInner"
                             , _tilex + wrapgap
                             , _tiley + wrapgap
                             , _tilez*tile_number + tilegap*(tile_number-1) + wrapgap
                             , 0, 0 );
  G4VSolid* wrapbox = new G4SubtractionSolid( "WrapBox"
                                            , wrapOuter, wrapInner );
  G4VSolid* wraphole = new G4Tubs( "WrapHole"
                                 , 0, _pcb_radius
                                 , _tilez
                                 , 0, 2*pi );


  const G4ThreeVector offset( 0, 0, 0.5*(_tilez*tile_number + tilegap*(tile_number-1)) + wrapgap + 0.5*wrapthickness );

  return new G4SubtractionSolid( "WrapSolid"
                               , wrapbox, wraphole
                               , NULL, offset );
}

G4VSolid*
LYSimDetectorConstruction::ConstructHollowWrapGapSolid() const
{
  G4VSolid* wrapgapOuter
    = ConstructTrapazoidSolid( "WrapGapOuter"
                             , _tilex + 2*wrapgap
                             , _tiley + 2*wrapgap
                             , _tilez*tile_number + tilegap*(tile_number-1) + 2*wrapgap
                             , 0, 0 );
  G4VSolid* wrapgapInner
    = ConstructTrapazoidSolid( "WrapGapInner"
                             , _tilex 
                             , _tiley
                             , _tilez*tile_number + tilegap*(tile_number-1)
                             , 0, 0 );
  G4VSolid* wrapgapbox = new G4SubtractionSolid( "WrapGapBox"
                                            , wrapgapOuter, wrapgapInner );
  G4VSolid* wrapgaphole = new G4Tubs( "WrapGapHole"
                                 , 0, _pcb_radius
                                 , _tilez
                                 , 0, 2*pi );

  const G4ThreeVector offset( 0, 0, 0.5*(_tilez*tile_number + tilegap*(tile_number-1)) + wrapgap + 0.5*wrapthickness );

  return new G4SubtractionSolid( "WrapGapSolid"
                               , wrapgapbox, wrapgaphole
                               , NULL, offset );
}
G4VSolid*
LYSimDetectorConstruction::ConstructUnEvenHollowWrapGapSolid() const
{
  G4VSolid* wrapgapOuter
    = ConstructTrapazoidSolid( "WrapGapOuter"
                             , _tilex + 2*0.1*mm
                             , _tiley + 2*0.1*mm
                             , _tilez*tile_number + tilegap*(tile_number-1) + 2*wrapgap
                             , 0, 0 );
  G4VSolid* wrapgapInner
    = ConstructTrapazoidSolid( "WrapGapInner"
                             , _tilex
                             , _tiley
                             , _tilez*tile_number + tilegap*(tile_number-1)
                             , 0, 0 );

  const G4ThreeVector unevenoffset( 0, 0, -1*wrapgap+0.1*mm );
  G4VSolid* wrapgapbox = new G4SubtractionSolid( "WrapGapBox"
                                            , wrapgapOuter, wrapgapInner, NULL, unevenoffset );
  G4VSolid* wrapgaphole = new G4Tubs( "WrapGapHole"
                                 , 0, _pcb_radius
                                 , _tilez
                                 , 0, 2*pi );

  const G4ThreeVector offset( 0, 0, 0.5*(_tilez*tile_number + tilegap*(tile_number-1)) + wrapgap + 0.5*wrapthickness );

  return new G4SubtractionSolid( "WrapGapSolid"
                               , wrapgapbox, wrapgaphole
                               , NULL, offset );
}

G4VSolid*
LYSimDetectorConstruction::ConstructSphereDimpleSolid() const
{
  G4Sphere* solidsphere = new G4Sphere( "DimpleSphere"
                                      , 0, GetDimpleSizeRadius()
                                      , 0.,  2.*pi
                                      , pi/2, pi );

  G4Box* solidsub = new G4Box( "DimpleSphereSubBox"
                             , GetDimpleSizeRadius(), GetDimpleSizeRadius()
                             , GetDimpleSizeRadius() - GetDimpleIndent() );

  return new G4SubtractionSolid( "Dimple", solidsphere, solidsub );
}

G4VSolid*
LYSimDetectorConstruction::ConstructParabolicDimpleSolid() const
{
  const G4double SemiZ = ( 447317/288000 )*mm + _dimple_indent;
  const G4double SemiX = ( 447317/120000 )*mm;

  return new G4Ellipsoid( "DimpleParabolic",
    SemiX,
    SemiX,
    SemiZ,
    -SemiZ,// must be < SemiZ
    -SemiZ+_dimple_indent );// must be >-SemiZ
}

G4VSolid*
LYSimDetectorConstruction::ConstructPyramidDimpleSolid() const
{
  return new G4Trd( "DimplePyramid",
    0.0*mm, _dimple_radius,
    0.0*mm, _dimple_radius,
    _dimple_indent*0.5 );
}

G4VSolid*
LYSimDetectorConstruction::ConstructArbiDimpleSolid() const
{
    //borrowed some functions that not used in arbitrary shape construction 
    const double rr  = GetDimpleRadius();
    const double ind  = GetDimpleIndent();
    const double dimp_s  = GetDimpleSA();

    //derived from sphere LCE
    TF1 *f=new TF1("f","(TMath::Log((3-[0]-x)/[1]))/[2]/TMath::Log(2.7182818)",0,60);
    f->SetParameter(0,rr);
    f->SetParameter(1,ind);
    f->SetParameter(2,dimp_s);

    TF1 *f_up=new TF1("f_up","3 - [0] - [1]*exp([2]*x)",0,60);
    f_up->SetParameter(0,rr);
    f_up->SetParameter(1,ind);
    f_up->SetParameter(2,dimp_s);

    double x_min = f_up->Eval(12);
    double x_max = 3;
    const int x_n=601;
    G4double zLB[x_n] ;
    G4double rinLB[x_n];
    G4double routLB[x_n];
    zLB[0] = 0;
    rinLB[0] = 0;
    routLB[0] = f->Eval(x_min) *mm;
    for(int i=1;i<=(x_n-1);i++){
        zLB[i] =  x_min*mm + ((i-1)*(x_max-x_min)*mm/(x_n-1));
        rinLB[i] = 0;
        if(f->Eval(zLB[i])>0)
          routLB[i] = f->Eval(zLB[i]) *mm;
        else
          routLB[i] = 0;
        zLB[i]*=-1;
    }

    return new G4Polycone("LittleBoySolid",0., 360.*deg, x_n, zLB, rinLB, routLB);
}

double
LYSimDetectorConstruction::GetDimpleSizeRadius() const
{
  return 0.5*( ( _dimple_radius*_dimple_radius )/_dimple_indent
               + _dimple_indent );
}

G4ThreeVector
LYSimDetectorConstruction::CalcDimpleOffset() const
{
  if( _dimple_type == PYRAMID ){
    return G4ThreeVector( 0, 0, 0.5*_tilez-0.5*_dimple_indent );
  } else if( _dimple_type == PARABOLIC ){
    return G4ThreeVector( 0, 0, 0.5*_tilez+( 447317/288000 )*mm );
  } else {
    return G4ThreeVector( 0, 0
                        , 0.5*_tilez+GetDimpleSizeRadius()-GetDimpleIndent() );
  }
}

G4ThreeVector
LYSimDetectorConstruction::CalcSiPMDimpleOffset() const
{
  if( _dimple_type == PYRAMID ){
    return G4ThreeVector( 0, 0, 0.5*_dimple_indent - 0.5*_sipm_z );
  } else if( _dimple_type == PARABOLIC ){
    return G4ThreeVector( 0, 0, 0 );// TO-DO!!
  } else {
    return G4ThreeVector( 0, 0,
      -GetDimpleSizeRadius() + GetDimpleIndent()- 0.5*_sipm_z );
  }
}

// Additional geometry factors
double
LYSimDetectorConstruction::WorldHalfX() const
{
  return _tilex * 1.5;
}

double
LYSimDetectorConstruction::WorldHalfY() const
{
  return _tiley *1.5;
}

double
LYSimDetectorConstruction::WorldHalfZ() const
{
  return _tilez * 100;
}

double
LYSimDetectorConstruction::LocalTileZ( const double x, const double y ) const
{
  const double dimple_r  = GetDimpleRadius();
  const double dimple_i  = GetDimpleIndent();
  const double dimple_s  = GetDimpleSA();

  const double beam_r    = sqrt( x*x+y*y );
  const double big_r     = GetDimpleSizeRadius();
  const double thickness = GetTileZ()-Gettgap()*(Gettn()-1);
  const int type         = GetDimpleType();

//**********************8
//add arbi shape 
//**********************8
/*
    //derived from sphere LCE
    TF1 *f_up=new TF1("f_up","[0]+[1]*exp([2]*x)",0,60);
    f_up->SetParameter(0,dimple_r);
    f_up->SetParameter(1,dimple_i);
    f_up->SetParameter(2,dimple_s);
    // radius from 0 to 12.000408
    if(beam_r<12){
      if( f_up->Eval(beam_r) >= 0) return  thickness - 3.0 + f_up->Eval(beam_r);
      else return 0;
    }else{
      return thickness;
    }

*/

  if( beam_r <= dimple_r ){
    if( type == SPHERICAL ){
      return (thickness - ( sqrt( big_r*big_r - beam_r*beam_r ) - ( big_r - dimple_i ) ));
    } else if( type == PARABOLIC ){//  # Parabolic
      return thickness - abs( 1.35 * beam_r*beam_r - dimple_i );
    } else if( type == PYRAMID ){
      if( fabs( y ) <= fabs( x ) ){
        return thickness - ( dimple_i - ( dimple_i / dimple_r ) * fabs( x ) );
      } else {
        return thickness - ( dimple_i - ( dimple_i / dimple_r ) * fabs( y ) );
      }
    }  else {
      return thickness;
    }
  } else {
      return thickness;
  }
}


void
LYSimDetectorConstruction::SetTileAbsMult( const double mult )
{
  _absmult = mult;
  if(material==0){Update_flat_AbsLength( fTileMaterial, _absmult );}
  else Update_EJ200_AbsLength( fTileMaterial, _absmult );
}

void
LYSimDetectorConstruction::SetSurfaceDimpleSA( const double mult )
{
  _dimplesa = mult;
  Update_Dimple_SigmaAlpha( fdimpleSurface, _dimplesa );
}
void
LYSimDetectorConstruction::SetWrapReflect( const double r )
{
  // Add entries into properties table
  _wrap_reflect = r;
  static const unsigned nentries = 43;
  static double phoE[nentries]   = {3.444005361           *eV,3.350924135           *eV,3.262741921           *eV,3.179081872           *eV,3.099604825           *eV,3.024004707           *eV,2.952004595           *eV,2.883353326           *eV,2.817822568           *eV,2.755204289           *eV,2.695308543           *eV,2.637961553           *eV,2.583004021           *eV,2.530289653           *eV,2.47968386           *eV,2.431062608           *eV,2.384311404           *eV,2.339324396           *eV,2.296003574           *eV,2.254258055           *eV,2.214003446           *eV,2.175161281           *eV,2.1376585           *eV,2.101427           *eV,2.066403217           *eV,2.032527754           *eV,1.999745048           *eV,1.968003063           *eV,1.937253016           *eV,1.907449123           *eV,1.878548379           *eV,1.850510343           *eV,1.823296956           *eV,1.796872362           *eV,1.771202757           *eV,1.746256239           *eV,1.722002681           *eV,1.698413603           *eV,1.675462068           *eV,1.653122573           *eV,1.631370961           *eV,1.610184325           *eV,1.589540936           *eV};
  double reflectivity[nentries]  = {0.1034           *_wrap_reflect,0.200833333           *_wrap_reflect,0.6738           *_wrap_reflect,0.868166667           *_wrap_reflect,0.915033333           *_wrap_reflect,0.927133333           *_wrap_reflect,0.933433333           *_wrap_reflect,0.942033333           *_wrap_reflect,0.9484           *_wrap_reflect,0.952633333           *_wrap_reflect,0.9566           *_wrap_reflect,0.9614           *_wrap_reflect,0.965166667           *_wrap_reflect,0.9679           *_wrap_reflect,0.969833333           *_wrap_reflect,0.970033333           *_wrap_reflect,0.969833333           *_wrap_reflect,0.969133333           *_wrap_reflect,0.966133333           *_wrap_reflect,0.963333333           *_wrap_reflect,0.965766667           *_wrap_reflect,0.968666667           *_wrap_reflect,0.9705           *_wrap_reflect,0.9716           *_wrap_reflect,0.971566667           *_wrap_reflect,0.972433333           *_wrap_reflect,0.972333333           *_wrap_reflect,0.970033333           *_wrap_reflect,0.947833333           *_wrap_reflect,0.936066667           *_wrap_reflect,0.948633333           *_wrap_reflect,0.949633333           *_wrap_reflect,0.9501           *_wrap_reflect,0.9545           *_wrap_reflect,0.955933333           *_wrap_reflect,0.957233333           *_wrap_reflect,0.961033333           *_wrap_reflect,0.962366667           *_wrap_reflect,0.963766667           *_wrap_reflect,0.9656           *_wrap_reflect,0.9662           *_wrap_reflect,0.9662           *_wrap_reflect,0.963933333           *_wrap_reflect};

  static const unsigned nentries_cons = 2;
  static double phoE_cons[nentries_cons]   = {1*eV, 8*eV};
  double reflectivity_cons[nentries_cons]  = {0.985*_wrap_reflect,0.985*_wrap_reflect};
  G4MaterialPropertiesTable* table = fESROpSurface->GetMaterialPropertiesTable();
  if( table ){
    table->RemoveProperty( "REFLECTIVITY" );
    table->AddProperty( "REFLECTIVITY", phoE_cons, reflectivity_cons, nentries_cons );
//    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  } else {
    table = new G4MaterialPropertiesTable();
//    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
    table->AddProperty( "REFLECTIVITY", phoE_cons, reflectivity_cons, nentries_cons );
    fESROpSurface->SetMaterialPropertiesTable( table );
  }

  static const unsigned nentriesty = 35;
 double phoEty[nentriesty]  = {3.444005361           *eV,3.350924135           *eV,3.262741921           *eV,3.179081872           *eV,3.099604825           *eV,3.024004707           *eV,2.952004595           *eV,2.883353326           *eV,2.817822568           *eV,2.755204289           *eV,2.695308543           *eV,2.637961553           *eV,2.583004021           *eV,2.530289653           *eV,2.47968386           *eV,2.431062608           *eV,2.384311404           *eV,2.339324396           *eV,2.296003574           *eV,2.254258055           *eV,2.214003446           *eV,2.175161281           *eV,2.1376585           *eV,2.101427           *eV,2.066403217           *eV,2.032527754           *eV,1.999745048           *eV,1.968003063           *eV,1.937253016           *eV,1.907449123           *eV,1.878548379           *eV,1.850510343           *eV,1.823296956           *eV,1.796872362           *eV,1.771202757           *eV};
// double reflectivityty[nentriesty]  = {0.823738509    *_wrap_reflect,0.842890705    *_wrap_reflect,0.878426966    *_wrap_reflect,0.906639428    *_wrap_reflect,0.913942799    *_wrap_reflect,0.933105209    *_wrap_reflect,0.952277835    *_wrap_reflect,0.96784474    *_wrap_reflect,0.98340143    *_wrap_reflect,0.988937692    *_wrap_reflect,0.989948927    *_wrap_reflect,0.992747702    *_wrap_reflect,0.991940756    *_wrap_reflect,0.993881512    *_wrap_reflect,0.99587334    *_wrap_reflect,0.996905005    *_wrap_reflect,0.997037794    *_wrap_reflect,0.999918284    *_wrap_reflect,0.997334014    *_wrap_reflect, 0.998396323    *_wrap_reflect,0.998518897    *_wrap_reflect,1    *_wrap_reflect,0.99691522    *_wrap_reflect,0.995209397    *_wrap_reflect,0.999846782    *_wrap_reflect,1    *_wrap_reflect };
 double reflectivityty[nentriesty]  = {0.7439           *_wrap_reflect,0.7532           *_wrap_reflect,0.7648           *_wrap_reflect,0.7681           *_wrap_reflect,0.77055           *_wrap_reflect,0.77115           *_wrap_reflect,0.7718           *_wrap_reflect,0.7727           *_wrap_reflect,0.7748           *_wrap_reflect,0.7759           *_wrap_reflect,0.77585           *_wrap_reflect,0.77555           *_wrap_reflect,0.7742           *_wrap_reflect,0.77385           *_wrap_reflect,0.77285           *_wrap_reflect,0.77195           *_wrap_reflect,0.7716           *_wrap_reflect,0.77085           *_wrap_reflect,0.7703           *_wrap_reflect,0.7702           *_wrap_reflect,0.77065           *_wrap_reflect,0.7703           *_wrap_reflect,0.7705           *_wrap_reflect,0.7713           *_wrap_reflect,0.7708           *_wrap_reflect,0.7704           *_wrap_reflect,0.77235           *_wrap_reflect,0.77265           *_wrap_reflect,0.77315           *_wrap_reflect,0.7732           *_wrap_reflect,0.7715           *_wrap_reflect,0.7705           *_wrap_reflect,0.7715           *_wrap_reflect,0.7714           *_wrap_reflect,0.77105           *_wrap_reflect,};


  G4MaterialPropertiesTable* tablet = fTyvekSurface->GetMaterialPropertiesTable();
  if( tablet ){
    tablet->RemoveProperty( "REFLECTIVITY" );
    tablet->AddProperty( "REFLECTIVITY", phoEty, reflectivityty, nentriesty );
  } else {
    tablet = new G4MaterialPropertiesTable();
    tablet->AddProperty( "REFLECTIVITY", phoEty, reflectivityty, nentriesty );
    fTyvekSurface->SetMaterialPropertiesTable( tablet );
  }
}

void
LYSimDetectorConstruction::SetPCBReflect( const double r )
{
  // Add entries into properties table
//1st measuement, with white paper
/*  const int nentries = 43;
  double phoE[nentries]   = {360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780 };
  double reflectivity[nentries] = {15.11, 15.4,      19.75,   27.935,   38.715,   52.025,   59.16,   62.25,   63.19,   63.44,   63.325,   63.18,   62.705,   62.295,   62.21,   62.035,   61.695,   61.26,   60.87,   60.505,   60.035,   59.485,   59.015,   58.535,   57.89,   57.295,   57.07,   56.525,   56.09,   55.68,   55.07,   54.65,   54.18,   53.72,   53.25,   52.83,   52.45,   52.195,   52.89,   54.16,   55.78,   57.59,   58.47};
*/
  const int nentries = 43;
  double phoE[nentries]   = {360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780};
//WSS
//  double reflectivity[nentries] = {13.66           ,26.35           ,43.59           ,51.13           ,57.46           ,61.25           ,63.26           ,63.96           ,63.98           ,63.74           ,63.39           ,62.83           ,62.01           ,61.26           ,60.56           ,59.84           ,59.13           ,58.43           ,57.77           ,57.15           ,56.52           ,55.85           ,55.25           ,54.69           ,54.14           ,53.43           ,52.98           ,52.35           ,51.93           ,51.89           ,52.72           ,53.52           ,54.19           ,54.76           ,55.26           ,55.71           ,56.21           ,56.76           ,57.54           ,58.08           ,58.46           ,58.64           ,58.83           };
//another measurment 
  double reflectivity[nentries] ={17.835           ,31.91           ,50.7           ,57.9           ,62.72           ,66.185           ,68.02           ,68.82           ,68.805           ,68.68           ,68.27           ,67.67           ,66.845           ,66.13           ,65.36           ,64.57           ,63.84           ,63.12           ,62.455           ,61.84           ,61.17           ,60.49           ,59.875           ,59.25           ,58.635           ,57.92           ,57.475           ,56.84           ,56.33           ,56.365           ,57.255           ,57.98           ,58.68           ,59.265           ,59.87           ,60.395           ,60.83           ,61.21           ,62.135           ,62.735           ,63.2           ,63.525           ,63.57           };
//blacktape
  double reflectivity2[nentries] ={5.4725           ,5.6975           ,5.15           ,5.085           ,5.095           ,5.19           ,5.1575           ,5.05           ,5.0275           ,5.0325           ,4.995           ,4.955           ,4.97           ,4.955           ,4.9475           ,4.945           ,4.915           ,4.8525           ,4.88           ,4.885           ,4.8875           ,4.925           ,4.885           ,4.86           ,4.845           ,4.93           ,4.755           ,4.665           ,4.68           ,4.6575           ,4.5125           ,4.6425           ,4.5375           ,4.47           ,4.42           ,4.3775           ,4.31           ,4.3475           ,4.3875           ,4.36           ,4.33           ,4.255           ,4.4525};

  _pcb_reflectivity = r;
  for(int i=0;i<nentries;i++){
    phoE[i] = 1239.84193*eV/phoE[i];
    reflectivity[i] = reflectivity[i]/100.0;
    reflectivity2[i] = reflectivity2[i]/100.0;
  }
  
  G4MaterialPropertiesTable* table = fPCBSurface->GetMaterialPropertiesTable();
  if( table ){
    table->RemoveProperty( "REFLECTIVITY" );
    if(_pcb_reflectivity<0.5) table->AddProperty( "REFLECTIVITY", phoE, reflectivity2, nentries );
    else table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  } else {
    table = new G4MaterialPropertiesTable();
    if(_pcb_reflectivity<0.5) table->AddProperty( "REFLECTIVITY", phoE, reflectivity2, nentries );
    else table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
    fPCBSurface->SetMaterialPropertiesTable( table );
  }


}




