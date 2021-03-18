#ifndef LYSimDetectorConstruction_h
#define LYSimDetectorConstruction_h 1

#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimPMTSD.hh"
#else
#include "LYSimPMTSD.hh"
#endif

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class LYSimDetectorMessenger;
class G4Box;
class G4LogicalVolume;
class G4Material;
class G4OpticalSurface;
class LYSimPMTSD;

class LYSimDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  LYSimDetectorConstruction();
  virtual ~LYSimDetectorConstruction();
  virtual G4VPhysicalVolume* Construct();
  void                       UpdateGeometry();

  // Get/Set functions of tile geometry parameters
  inline void
  SetTileZ ( const G4double x ){ _tilez = x;}
  inline G4double
  GetTileZ() const {return _tilez;}
  inline void
  SetTileX ( const G4double x ){_tilex = x;}
  inline G4double
  GetTileX() const {return _tilex;}
  inline void
  SetTileY ( const G4double y ){_tiley = y;}
  inline G4double
  GetTileY() const {return _tiley;}

  // Get/Set functions of dimple geometry parameters
  
  inline double
  GetDimpleRadius() const { return _dimple_radius; }
  inline void
  SetDimpleRadius( const double x ){ _dimple_radius = x; }
  inline double
  GetDimpleIndent() const { return _dimple_indent; }
  inline void
  SetDimpleIndent( const double x ){ _dimple_indent = x; }
  int
  GetDimpleType() const {return _dimple_type; }
  inline void
  SetDimpleType( const int x ){ _dimple_type = x; }
  

  // Get/Set functions for SiPM geometry parameters
  inline double
  GetSiPMX() const {return _sipm_x; }
  inline void
  SetSiPMX( const double x ){ _sipm_x = x; }
  inline double
  GetSiPMY() const { return _sipm_y; }
  inline void
  SetSiPMY( const double x ){ _sipm_y = x; }
  inline double
  GetSiPMRim() const {return _sipm_rimwidth; }
  inline void
  SetSiPMRim( const double x ){ _sipm_rimwidth = x; }
  inline double
  GetSiPMStand() const { return _sipm_standz; }
  inline void
  SetSiPMStand( const double x ){ _sipm_standz = x; }

  enum  DIMPLE_TYPE
  {
    SPHERICAL,
    ELLIPSOID,
    FLAT_DOME,
    CYLINDRICAL
  };

  // Additional geometric factors.
  double WorldHalfX() const;
  double WorldHalfY() const;
  double WorldHalfZ() const;

  double LocalTileZ( const double x, const double y ) const;

  // Material updating functions
  void SetTileAbsMult( const double x );
  void SetTileScintillation( const double x );
  inline double
  GetTileAbsMult() const { return _absmult; }
  double
  GetTileScintillation() const { return _ScintiN; }

  void SetWrapReflect( const double x );
  inline double
  GetWrapReflect() const { return _wrap_reflect; }

  void SetTileAlpha( const double x );
  inline double
  GetTileAlpha(){ return _tile_alpha; }

  void SetDimpleAlpha( const double x );
  inline double
  GetDimpleAlpha(){ return _dimple_alpha; }

  void
  SetPCBRadius( const double x ){ _pcb_radius = x; }
  inline double
  GetPCBRadius() const { return _pcb_radius; }

  void SetPCBReflect( const double x );
  inline double
  GetPCBReflect() const { return _pcb_reflectivity; }

  // Accesing Custom Materials
  inline G4Material*
  GetEJ200(){ return fEJ200; }
  inline G4OpticalSurface*
  GetSiPMSurface(){ return fSiPMSurface; }


  //wls
  void
  Set_handwrap( const bool x ){ _handwrap = x; }
  inline bool
  Get_handwrap() const { return _handwrap; }
  void
  SetHoleRadius( const double x ){ _hole_radius = x; }
  inline double
  GetHoleRadius() const { return _hole_radius; }
  void
  SetHoleX1( const double x ){ _hole_x1 = x; }
  inline double
  GetHoleX1() const { return _hole_x1; }
  void
  SetHoleX2( const double x ){ _hole_x2 = x; }
  inline double
  GetHoleX2() const { return _hole_x2; }
  void
  SetFiberR( const double x ){ _WLSfiberR = x; }
  inline double
  GetFiberR() const { return _WLSfiberR; }
  void
  SetFiber_clad_thick( const double x ){ _WLSfiber_clad_thick = x; }
  inline double
  GetFiber_clad_thick() const { return _WLSfiber_clad_thick; }
  void
  SetFiberZ( const double x ){ _WLSfiberZ = x*CLHEP::m; }
  inline double
  GetFiberZ() const { return _WLSfiberZ; }
  void
  SetFiberZoff( const double x ){ _WLS_zoff = x*CLHEP::m; }
  inline double
  GetFiberZoff() const { return _WLS_zoff; }


private:
  // Subfunctions for cleaner code
  G4VSolid* ConstructTrapazoidSolid(
    const G4String& name,
    double          x,// Maximum Length along z direction.
    double          y,// Length along y direction
    double          z,// length along z direction
    double          indent_x1,// Indent in x direction
    double          indent_x2// indent of corner
    ) const;

  G4VSolid* ConstructHollowWrapSolid() const;
  G4VSolid* ConstructHollowWrapCladSolid() const;

  G4VSolid* ConstructDimpleSubtract() const;

  G4VSolid* ConstructSphereDimpleSolid() const;
  G4VSolid* ConstructFlatDomeDimpleSolid() const;
  G4VSolid* ConstructCylindricalDimpleSolid() const;
  G4VSolid* ConstructEllipsoidDimpleSolid() const;

  G4ThreeVector CalcDimpleOffset() const;
  double        LocalTileZSpherical( const double r ) const;
  double        LocalTileZFlatDome( const double r ) const;
  double        LocalTileZCylindrical( const double r ) const;
  double        LocalTileZEllipsoid( const double r ) const;

  double GetDimpleSphereRadius() const;

  // Pointer to detector messenger class
  LYSimDetectorMessenger* fdetectorMessenger;

  // Pointers to materials
  G4Material* fAir;
  G4Material* fEJ200;
  G4Material* fBialkali;
  G4Material* fEpoxy;
  G4Material* fResin;

  // Pointers to surfaces
  G4OpticalSurface* fESROpSurface;
  G4OpticalSurface* fIdealPolishedOpSurface;
  G4OpticalSurface* fTileBulkSurface;
  G4OpticalSurface* fTileDimpleSurface;
  G4OpticalSurface* fIdealWhiteOpSurface;
  G4OpticalSurface* fSiPMSurface;
  G4OpticalSurface* fPCBSurface;

  // Pointers for access to Sensitive Detector
  static LYSimPMTSD* fPMTSD;

  // Tile Geometric definitions
  double _tilex;
  double _tiley;
  double _tilez;
  double _tile_x1;
  double _tile_x2;

  // Wrapping Geometry definitions
  double wrapgap;// distance between tile and wrapping
  double wrapthickness;

  // Absorption length multiplier
  double _absmult;
  double _wrap_reflect;
  double _tile_alpha;
  double _dimple_alpha;

  // Dimple variables;
  double _dimple_indent;// Depth of dimple
  double _dimple_radius;// Dimple radius
  int _dimple_type;

  // Parameters for plate behind SiPM
  double _pcb_radius;
  double _pcb_reflectivity;

  double _sipm_deadwidth;
  double _sipm_x;
  double _sipm_y;
  double _sipm_z;
  double _sipm_standz;
  double _sipm_rimwidth;
  double _sipm_glasswidth;

  //wls
  bool _handwrap;
  double _hole_radius;
  double _hole_x1;
  double _hole_x2;
  double _WLSfiberR;
  double _WLSfiberZ;
  double _WLSfiber_clad_thick;
  double _WLS_zoff;
  double _ScintiN;
  G4Material* mfiber;
  G4Material* mfiber_clad;
  G4Material* fcoating;
  G4Material* fholemat;
  G4OpticalSurface* fTiO2Surface;
  G4OpticalSurface* opSurface;
};

#endif
