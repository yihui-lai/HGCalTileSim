#include "G4NistManager.hh"
#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/Materials.hh"
#else
#include "Materials.hh"
#endif

using namespace CLHEP;

// static list of common materials (cannot be const)
static G4Element* C  = G4NistManager::Instance()->FindOrBuildElement( "C" );
static G4Element* H  = G4NistManager::Instance()->FindOrBuildElement( "H" );
static G4Element* Si = G4NistManager::Instance()->FindOrBuildElement( "Si" );
static G4Element* O  = G4NistManager::Instance()->FindOrBuildElement( "O" );
static G4Element* K  = G4NistManager::Instance()->FindOrBuildElement( "K" );
static G4Element* Sb = G4NistManager::Instance()->FindOrBuildElement( "Sb" );
static G4Element* Rb = G4NistManager::Instance()->FindOrBuildElement( "Rb" );
static G4Element* Cs = G4NistManager::Instance()->FindOrBuildElement( "Cs" );
static G4Element* Lu = G4NistManager::Instance()->FindOrBuildElement( "Lu" );
static G4Element* Y  = G4NistManager::Instance()->FindOrBuildElement( "Y" );
static G4Element* Ce = G4NistManager::Instance()->FindOrBuildElement( "Ce" );
static G4Element* La = G4NistManager::Instance()->FindOrBuildElement( "La" );
static G4Element* Br = G4NistManager::Instance()->FindOrBuildElement( "Br" );
static G4Element* Na = G4NistManager::Instance()->FindOrBuildElement( "Na" );
static G4Element* I  = G4NistManager::Instance()->FindOrBuildElement( "I" );
static G4Element* Tl = G4NistManager::Instance()->FindOrBuildElement( "Tl" );
static G4Element* Al = G4NistManager::Instance()->FindOrBuildElement( "Al" );

static const unsigned flatentries = 2;
static const double minenergy     = 1.0*eV;
static const double maxenergy     = 8.0*eV;

// Bialkali definition:
// Ref: http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2013/03/11/12.39-85121-chDetectorConstruction.cc
G4Material*
Make_Bialkali()
{
  G4Material* material
    = new G4Material( "Bialkali", 4.28*g/cm3, 3, kStateSolid );
  material->AddElement( K,  13.3*perCent );
  material->AddElement( Cs, 45.2*perCent );
  material->AddElement( Sb, 41.5*perCent );
  return material;
}

G4Material*
Make_Epoxy()
{
  G4Material* epoxy = new G4Material( "Epoxy", 3.75*g/cm3, 2, kStateSolid );
  epoxy->AddElement( Al, 40*perCent );
  epoxy->AddElement( O,  60*perCent );
  double photonE[flatentries]      = {minenergy, maxenergy};
  double refrac_idx[flatentries]   = {1.50, 1.50};
  double Reflectivity[flatentries] = {0.0, 0.0};
  // Will obe overwritten in surface property defintiion
  double abs_length[flatentries] = {10*m, 10*m};

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();
  table->AddProperty( "RINDEX",  photonE, refrac_idx,   flatentries );
  table->AddProperty( "REFLECTIVITY", photonE, Reflectivity, flatentries );
  table->AddProperty( "ABSLENGTH",    photonE, abs_length,   flatentries );

  epoxy->SetMaterialPropertiesTable( table );
  return epoxy;
}

G4Material*
Make_Custom_Air()
{
  static double photonE[flatentries]    = {minenergy, maxenergy};
  static double refrac_idx[flatentries] = {1.0, 1.0};

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();
  table->AddProperty( "RINDEX", photonE, refrac_idx, flatentries );

  G4Material* material
    = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );
  material->SetMaterialPropertiesTable( table );
  return material;
}

G4Material*
Make_Resin()
{
  static double photonE[flatentries]    = {minenergy, maxenergy};
  static double refrac_idx[flatentries] = {1.52, 1.52};
  static double abs_length[flatentries] = {10*m, 10*m};

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();
  table->AddProperty( "RINDEX",    photonE, refrac_idx, flatentries );
  table->AddProperty( "ABSLENGTH", photonE, abs_length, flatentries );

  G4Material* material
    = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Pyrex_Glass" );
  material->SetMaterialPropertiesTable( table );

  return material;
}

G4Material*
Make_EJ200()
{
  G4Material* material = new G4Material( "EJ200", 1.023*g/cm3, 2, kStateSolid );
  material->AddElement( C, 91.53*perCent );
  material->AddElement( H, 8.47*perCent );

  static double photonE[flatentries]    = {minenergy, maxenergy};
  static double refrac_idx[flatentries] = { 1.58, 1.58 };

  const G4int NUMENTRIES_1 = 139;
  static double FAST_Energy[NUMENTRIES_1] = {3.13 * eV, 3.12 * eV, 3.11 * eV, 3.10 * eV, 3.09 * eV, 3.07 * eV, 3.07 * eV, 3.06 * eV, 3.05 * eV, 3.05 * eV, 3.05 * eV, 3.04 * eV, 3.04 * eV,    3.03 * eV, 3.03 * eV, 3.02 * eV, 3.02 * eV, 3.02 * eV, 3.01 * eV, 3.01 * eV, 3.00 * eV, 3.00 * eV, 3.00 * eV,    3.00 * eV, 2.99 * eV, 2.99 * eV, 2.99 * eV, 2.98 * eV, 2.98 * eV, 2.97 * eV, 2.97 * eV, 2.98 * eV, 2.98 * eV,    2.97 * eV, 2.97 * eV, 2.96 * eV, 2.96 * eV, 2.96 * eV, 2.95 * eV, 2.95 * eV, 2.95 * eV, 2.94 * eV, 2.94 * eV,    2.93 * eV, 2.93 * eV, 2.92 * eV, 2.91 * eV, 2.90 * eV, 2.89 * eV, 2.89 * eV, 2.88 * eV, 2.88 * eV, 2.87 * eV,    2.86 * eV, 2.85 * eV, 2.85 * eV, 2.84 * eV, 2.84 * eV, 2.83 * eV, 2.83 * eV, 2.82 * eV, 2.82 * eV, 2.81 * eV,    2.81 * eV, 2.80 * eV, 2.80 * eV, 2.79 * eV, 2.78 * eV, 2.78 * eV, 2.77 * eV, 2.77 * eV, 2.76 * eV, 2.76 * eV,    2.75 * eV, 2.74 * eV, 2.74 * eV, 2.73 * eV, 2.72 * eV, 2.72 * eV, 2.71 * eV, 2.70 * eV, 2.70 * eV, 2.69 * eV,    2.68 * eV, 2.67 * eV, 2.67 * eV, 2.67 * eV, 2.66 * eV, 2.66 * eV, 2.65 * eV, 2.65 * eV, 2.64 * eV, 2.64 * eV,    2.63 * eV, 2.63 * eV, 2.62 * eV, 2.62 * eV, 2.61 * eV, 2.60 * eV, 2.60 * eV, 2.59 * eV, 2.59 * eV, 2.58 * eV,    2.58 * eV, 2.57 * eV, 2.56 * eV, 2.55 * eV, 2.55 * eV, 2.54 * eV, 2.53 * eV, 2.52 * eV, 2.52 * eV, 2.51 * eV,    2.51 * eV, 2.50 * eV, 2.49 * eV, 2.49 * eV, 2.48 * eV, 2.48 * eV, 2.46 * eV, 2.43 * eV, 2.41 * eV, 2.38 * eV,    2.36 * eV, 2.34 * eV, 2.32 * eV, 2.30 * eV, 2.28 * eV, 2.25 * eV, 2.23 * eV, 2.21 * eV, 2.19 * eV, 2.18 * eV,    2.16 * eV, 2.14 * eV, 2.12 * eV, 2.10 * eV, 2.07 * eV, 1.91 * eV};
  static double FAST_COMPONENT[NUMENTRIES_1] = {0, 0.0008, 0.0065, 0.0089, 0.0185, 0.0355, 0.0573, 0.0839, 0.1081, 0.1274, 0.1444, 0.171, 0.1927, 0.2218,       0.2508, 0.2677, 0.2968, 0.3185, 0.35, 0.379, 0.4056, 0.4323, 0.4565, 0.4903,       0.5218, 0.5532, 0.5798, 0.6065, 0.6282, 0.7177, 0.6984, 0.679, 0.6548, 0.7419,       0.7637, 0.7927, 0.8218, 0.8483, 0.8774, 0.9065, 0.9355, 0.9597, 0.979, 0.9887,       0.9983, 1, 0.996, 0.9887, 0.979, 0.9645, 0.95, 0.9355, 0.9161, 0.8968, 0.8774,       0.8605, 0.8387, 0.8145, 0.7903, 0.7734, 0.7565, 0.7371, 0.7177, 0.6935, 0.6766,       0.6597, 0.6452, 0.6258, 0.6089, 0.5919, 0.5774, 0.5629, 0.5508, 0.5315, 0.5145,       0.5, 0.4855, 0.471, 0.4565, 0.4468, 0.4371, 0.4274, 0.4177, 0.4081, 0.396, 0.379,       0.3887, 0.3669, 0.3548, 0.3403, 0.3258, 0.3016, 0.3113, 0.2823, 0.2629, 0.2484,       0.2339, 0.2242, 0.2097, 0.2, 0.1855, 0.1758, 0.1661, 0.154, 0.1468, 0.1323,       0.1226, 0.1129, 0.1032, 0.096, 0.0887, 0.0839, 0.0766, 0.0742, 0.0669, 0.0597,       0.0573, 0.0548, 0.0536, 0.0402, 0.0302, 0.0226, 0.017, 0.0127, 0.0096, 0.0072,       0.0054, 0.004, 0.003, 0.0023, 0.0017, 0.0013, 0.001, 0.0007, 0.0005, 0.0004, 0.0003, 0.0002, 0};

  // Add entries into properties table
  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("RINDEX", photonE, refrac_idx, flatentries);
  table->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);
  table->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
  table->AddConstProperty("RESOLUTIONSCALE", 1.0);
  table->AddConstProperty("FASTTIMECONSTANT", 2.5 * ns);
  table->AddConstProperty("YIELDRATIO", 1.0);
  table->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.5 * ns);

  material->SetMaterialPropertiesTable( table );

  // FIXME: Set the Birks Constant for the EJ200 scintillator
  material->GetIonisation()->SetBirksConstant( 0.126*mm/MeV );

  Update_EJ200_AbsLength( material, 1 );
  Update_EJ200_Scinti( material, 10 );
  return material;
}
void Update_Y11_timeconstant(G4Material* material, const double mult)
{
  G4MaterialPropertiesTable *tabler = material->GetMaterialPropertiesTable();
  tabler->RemoveConstProperty("WLSTIMECONSTANT");
  tabler->AddConstProperty("WLSTIMECONSTANT", mult * ns); //accounts for any time delay which may occur between absorption and re-emission of the photon, defalt delta
}

void
Update_EJ200_Scinti(G4Material* material, const double mult)
{

  G4MaterialPropertiesTable *tabler = material->GetMaterialPropertiesTable();
  tabler->RemoveConstProperty("SCINTILLATIONYIELD");
  tabler->AddConstProperty("SCINTILLATIONYIELD", mult / MeV);

}

void
Update_EJ200_AbsLength( G4Material* material, const double mult )
{
  //umd measurement
  
  //static const unsigned nentries = 401;
  //static double photonE[nentries] = {2.066403217 * eV, 2.069852972 * eV, 2.073314264 * eV, 2.076787152 * eV, 2.080271695 * eV, 2.08376795 * eV, 2.087275976 * eV, 2.090795835 * eV, 2.094327584 * eV, 2.097871286 * eV, 2.101427 * eV, 2.104994788 * eV, 2.108574711 * eV, 2.112166831 * eV, 2.115771212 * eV, 2.119387915 * eV, 2.123017003 * eV, 2.126658542 * eV, 2.130312595 * eV, 2.133979225 * eV, 2.1376585 * eV, 2.141350484 * eV, 2.145055242 * eV, 2.148772842 * eV, 2.152503351 * eV, 2.156246835 * eV, 2.160003362 * eV, 2.163773002 * eV, 2.167555822 * eV, 2.171351891 * eV, 2.175161281 * eV, 2.17898406 * eV, 2.182820299 * eV, 2.186670071 * eV, 2.190533445 * eV, 2.194410496 * eV, 2.198301294 * eV, 2.202205915 * eV, 2.206124431 * eV, 2.210056916 * eV, 2.214003446 * eV, 2.217964097 * eV, 2.221938943 * eV, 2.225928061 * eV, 2.229931529 * eV, 2.233949423 * eV, 2.237981823 * eV, 2.242028807 * eV, 2.246090453 * eV, 2.250166842 * eV, 2.254258055 * eV, 2.258364171 * eV, 2.262485274 * eV, 2.266621444 * eV, 2.270772766 * eV, 2.274939321 * eV, 2.279121195 * eV, 2.283318471 * eV, 2.287531236 * eV, 2.291759575 * eV, 2.296003574 * eV, 2.300263321 * eV, 2.304538903 * eV, 2.30883041 * eV, 2.313137929 * eV, 2.317461551 * eV, 2.321801367 * eV, 2.326157467 * eV, 2.330529944 * eV, 2.334918889 * eV, 2.339324396 * eV, 2.34374656 * eV, 2.348185473 * eV, 2.352641233 * eV, 2.357113935 * eV, 2.361603676 * eV, 2.366110553 * eV, 2.370634665 * eV, 2.375176111 * eV, 2.37973499 * eV, 2.384311404 * eV, 2.388905453 * eV, 2.393517239 * eV, 2.398146867 * eV, 2.402794438 * eV, 2.407460058 * eV, 2.412143833 * eV, 2.416845867 * eV, 2.42156627 * eV, 2.426305147 * eV, 2.431062608 * eV, 2.435838762 * eV, 2.44063372 * eV, 2.445447594 * eV, 2.450280494 * eV, 2.455132535 * eV, 2.460003829 * eV, 2.464894493 * eV, 2.469804641 * eV, 2.474734391 * eV, 2.47968386 * eV, 2.484653166 * eV, 2.48964243 * eV, 2.494651771 * eV, 2.49968131 * eV, 2.504731172 * eV, 2.509801478 * eV, 2.514892353 * eV, 2.520003923 * eV, 2.525136314 * eV, 2.530289653 * eV, 2.53546407 * eV, 2.540659693 * eV, 2.545876653 * eV, 2.551115082 * eV, 2.556375113 * eV, 2.56165688 * eV, 2.566960518 * eV, 2.572286162 * eV, 2.57763395 * eV, 2.583004021 * eV, 2.588396514 * eV, 2.593811569 * eV, 2.599249329 * eV, 2.604709937 * eV, 2.610193537 * eV, 2.615700274 * eV, 2.621230296 * eV, 2.62678375 * eV, 2.632360786 * eV, 2.637961553 * eV, 2.643586205 * eV, 2.649234893 * eV, 2.654907773 * eV, 2.660605 * eV, 2.666326731 * eV, 2.672073125 * eV, 2.677844341 * eV, 2.683640541 * eV, 2.689461887 * eV, 2.695308543 * eV, 2.701180675 * eV, 2.70707845 * eV, 2.713002035 * eV, 2.718951601 * eV, 2.724927319 * eV, 2.730929361 * eV, 2.736957903 * eV, 2.743013119 * eV, 2.749095188 * eV, 2.755204289 * eV, 2.761340601 * eV, 2.767504308 * eV, 2.773695593 * eV, 2.779914641 * eV, 2.78616164 * eV, 2.792436779 * eV, 2.798740248 * eV, 2.80507224 * eV, 2.811432948 * eV, 2.817822568 * eV, 2.824241298 * eV, 2.830689338 * eV, 2.837166888 * eV, 2.843674151 * eV, 2.850211333 * eV, 2.856778641 * eV, 2.863376282 * eV, 2.870004468 * eV, 2.876663411 * eV, 2.883353326 * eV, 2.890074429 * eV, 2.896826939 * eV, 2.903611077 * eV, 2.910427066 * eV, 2.917275129 * eV, 2.924155495 * eV, 2.931068392 * eV, 2.938014052 * eV, 2.944992708 * eV, 2.952004595 * eV, 2.959049952 * eV, 2.966129019 * eV, 2.973242038 * eV, 2.980389255 * eV, 2.987570916 * eV, 2.994787271 * eV, 3.002038571 * eV, 3.009325073 * eV, 3.016647032 * eV, 3.024004707 * eV, 3.031398362 * eV, 3.03882826 * eV, 3.046294668 * eV, 3.053797857 * eV, 3.061338099 * eV, 3.068915668 * eV, 3.076530844 * eV, 3.084183905 * eV, 3.091875137 * eV, 3.099604825 * eV, 3.107373258 * eV, 3.115180729 * eV, 3.123027531 * eV, 3.130913965 * eV, 3.138840329 * eV, 3.146806929 * eV, 3.154814071 * eV, 3.162862066 * eV, 3.170951228 * eV, 3.179081872 * eV, 3.187254319 * eV, 3.195468892 * eV, 3.203725917 * eV, 3.212025725 * eV, 3.220368649 * eV, 3.228755026 * eV, 3.237185196 * eV, 3.245659503 * eV, 3.254178294 * eV, 3.262741921 * eV, 3.271350739 * eV, 3.280005106 * eV, 3.288705385 * eV, 3.297451941 * eV, 3.306245147 * eV, 3.315085374 * eV, 3.323973003 * eV, 3.332908414 * eV, 3.341891995 * eV, 3.350924135 * eV, 3.36000523 * eV, 3.369135679 * eV, 3.378315886 * eV, 3.387546257 * eV, 3.396827205 * eV, 3.406159148 * eV, 3.415542507 * eV, 3.424977707 * eV, 3.43446518 * eV, 3.444005361 * eV, 3.453598691 * eV, 3.463245615 * eV, 3.472946583 * eV, 3.482702051 * eV, 3.492512479 * eV, 3.502378333 * eV, 3.512300085 * eV, 3.52227821 * eV, 3.532313191 * eV, 3.542405514 * eV, 3.552555673 * eV, 3.562764167 * eV, 3.573031499 * eV, 3.583358179 * eV, 3.593744725 * eV, 3.604191657 * eV, 3.614699504 * eV, 3.625268801 * eV, 3.635900088 * eV, 3.646593912 * eV, 3.657350826 * eV, 3.668171391 * eV, 3.679056172 * eV, 3.690005744 * eV, 3.701020687 * eV, 3.712101587 * eV, 3.723249039 * eV, 3.734463645 * eV, 3.745746012 * eV, 3.757096758 * eV, 3.768516505 * eV, 3.780005884 * eV, 3.791565535 * eV, 3.803196104 * eV, 3.814898246 * eV, 3.826672623 * eV, 3.838519907 * eV, 3.850440776 * eV, 3.862435919 * eV, 3.874506031 * eV, 3.886651818 * eV, 3.898873994 * eV, 3.911173281 * eV, 3.923550411 * eV, 3.936006127 * eV, 3.948541178 * eV, 3.961156326 * eV, 3.97385234 * eV, 3.98663 * eV, 3.999490097 * eV, 4.01243343 * eV, 4.025460812 * eV, 4.038573062 * eV, 4.051771013 * eV, 4.065055508 * eV, 4.078427401 * eV, 4.091887558 * eV, 4.105436854 * eV, 4.119076179 * eV, 4.132806433 * eV, 4.146628528 * eV, 4.160543389 * eV, 4.174551953 * eV, 4.188655169 * eV, 4.202854 * eV, 4.217149422 * eV, 4.231542423 * eV, 4.246034007 * eV, 4.260625189 * eV, 4.275317 * eV, 4.290110484 * eV, 4.305006701 * eV, 4.320006725 * eV, 4.335111643 * eV, 4.350322561 * eV, 4.365640599 * eV, 4.38106689 * eV, 4.396602589 * eV, 4.412248861 * eV, 4.428006893 * eV, 4.443877885 * eV, 4.459863058 * eV, 4.475963646 * eV, 4.492180906 * eV, 4.508516109 * eV, 4.524970547 * eV, 4.541545531 * eV, 4.55824239 * eV, 4.575062472 * eV, 4.592007148 * eV, 4.609077807 * eV, 4.626275858 * eV, 4.643602734 * eV, 4.661059887 * eV, 4.678648792 * eV, 4.696370947 * eV, 4.714227871 * eV, 4.732221107 * eV, 4.750352222 * eV, 4.768622808 * eV, 4.787034479 * eV, 4.805588876 * eV, 4.824287665 * eV, 4.843132539 * eV, 4.862125216 * eV, 4.881267441 * eV, 4.900560988 * eV, 4.920007659 * eV, 4.939609283 * eV, 4.95936772 * eV, 4.979284859 * eV, 4.999362621 * eV, 5.019602955 * eV, 5.040007846 * eV, 5.060579306 * eV, 5.081319385 * eV, 5.102230165 * eV, 5.12331376 * eV, 5.144572324 * eV, 5.166008042 * eV, 5.187623138 * eV, 5.209419874 * eV, 5.231400549 * eV, 5.2535675 * eV, 5.275923106 * eV, 5.298469786 * eV, 5.32121 * eV, 5.34414625 * eV, 5.367281082 * eV, 5.390617087 * eV, 5.4141569 * eV, 5.437903202 * eV, 5.461858722 * eV, 5.486026239 * eV, 5.510408578 * eV, 5.535008616 * eV, 5.559829283 * eV, 5.584873559 * eV, 5.61014448 * eV, 5.635645136 * eV, 5.661378676 * eV, 5.687348303 * eV, 5.713557281 * eV, 5.740008935 * eV, 5.766706651 * eV, 5.793653879 * eV, 5.820854131 * eV, 5.848310991 * eV, 5.876028104 * eV, 5.90400919 * eV, 5.932258038 * eV, 5.96077851 * eV, 5.989574541 * eV, 6.018650146 * eV, 6.048009415 * eV, 6.07765652 * eV, 6.107595714 * eV, 6.137831337 * eV, 6.168367811 * eV, 6.19920965 * eV};
  //static double Normabs_length[nentries] = {31.181,   31.096,   30.567,   30.233,   30.291,   29.862,   29.527,   29.45,   29.211,   29.061,   28.69,   28.088,   28.294,   28.091,   27.618,   28.104,   27.542,   27.323,   27.129,   27.653,   26.973,   26.923,   26.67,   26.931,   26.548,   26.396,   26.403,   26.134,   25.851,   25.951,   26.089,   25.929,   25.564,   25.623,   25.563,   25.485,   25.687,   25.809,   25.268,   25.195,   24.988,   24.614,   24.644,   24.803,   24.566,   24.498,   24.411,   24.271,   24.023,   24.098,   24.122,   23.971,   23.884,   23.866,   23.49,   23.509,   23.405,   23.49,   23.351,   23.137,   23.402,   22.929,   22.814,   22.715,   22.599,   22.536,   22.572,   22.518,   22.326,   22.319,   22.32,   22.021,   21.932,   21.927,   21.791,   21.634,   21.455,   21.51,   21.438,   21.261,   21.303,   21.281,   21.15,   21.187,   21.022,   20.74,   20.729,   20.691,   20.304,   20.552,   20.285,   20.231,   20.163,   20.051,   19.994,   19.79,   19.774,   19.513,   19.392,   19.338,   19.417,   19.355,   19.041,   19.047,   18.909,   18.735,   18.812,   18.942,   18.832,   18.674,   18.641,   18.5,   18.388,   18.326,   18.173,   18.227,   18.094,   17.866,   17.785,   17.564,   17.404,   17.232,   17.14,   17.001,   17.06,   16.868,   16.789,   16.739,   16.493,   16.506,   16.455,   16.387,   16.189,   16.08,   15.97,   15.977,   15.915,   15.75,   15.768,   15.564,   15.627,   15.467,   15.371,   15.295,   15.134,   15.058,   15.017,   14.868,   14.888,   14.682,   14.731,   14.59,   14.547,   14.439,   14.272,   14.243,   14.185,   14.058,   13.956,   13.885,   13.71,   13.649,   13.49,   13.616,   13.46,   13.339,   13.269,   13.174,   12.975,   12.98,   12.814,   12.783,   12.629,   12.353,   12.268,   12.09,   11.791,   11.503,   11.258,   10.91,   10.509,   10.089,   9.6212,   9.0457,   8.4357,   7.6776,   6.8687,   5.9825,   5.1017,   4.2291,   3.4321,   2.7202,   2.1268,   1.6411,   1.2669,   0.97412,   0.7556,   0.59057,   0.46697,   0.37392,   0.30354,   0.25051,   0.21089,   0.18138,   0.16275,   0.15173,   0.14624,   0.14404,   0.14284,   0.14264,   0.14216,   0.14216,   0.14199,   0.14173,   0.14185,   0.14195,   0.14174,   0.1419,   0.14166,   0.14232,   0.14134,   0.14213,   0.14161,   0.14201,   0.14184,   0.14191,   0.14196,   0.14197,   0.14154,   0.142,   0.14182,   0.14237,   0.14242,   0.14244,   0.14232,   0.14139,   0.14258,   0.14335,   0.14246,   0.14281,   0.14224,   0.14272,   0.14267,   0.14283,   0.14293,   0.14267,   0.14331,   0.14136,   0.14267,   0.14386,   0.14317,   0.14181,   0.14376,   0.14487,   0.14246,   0.14306,   0.1433,   0.14296,   0.1428,   0.14321,   0.14337,   0.14395,   0.14233,   0.14372,   0.14254,   0.14247,   0.14177,   0.14169,   0.14258,   0.14229,   0.141,   0.13997,   0.14017,   0.13996,   0.1382,   0.14058,   0.14048,   0.13797,   0.13932,   0.13866,   0.13838,   0.13915,   0.1371,   0.1388,   0.13793,   0.13734,   0.13692,   0.13784,   0.13922,   0.13762,   0.13928,   0.13739,   0.13833,   0.13861,   0.13848,   0.13891,   0.1391,   0.13859,   0.13897,   0.13855,   0.13906,   0.13819,   0.1386,   0.13775,   0.13894,   0.13863,   0.13822,   0.13831,   0.13826,   0.13879,   0.13858,   0.13858,   0.13793,   0.13843,   0.13773,   0.13776,   0.13804,   0.13793,   0.13772,   0.13729,   0.13829,   0.13699,   0.13661,   0.13689,   0.13648,   0.13596,   0.13628,   0.13571,   0.1363,   0.13601,   0.13647,   0.13546,   0.13538,   0.13565,   0.1354,   0.13567,   0.13611,   0.13546,   0.13582,   0.13621,   0.13637,   0.13596,   0.13636,   0.13683,   0.13659,   0.13631,   0.13629,   0.13629,   0.13618,   0.13639,   0.13599,   0.13619,   0.13559,   0.13579,   0.13556,   0.13566,   0.13524,   0.13558,   0.13515,   0.13571,   0.13521,   0.13513,   0.13539,   0.13499,   0.13505,   0.13467,   0.13452,   0.13418,   0.13386,   0.1327,   0.13145,   0.12959,   0.12899,   0.127,   0.12591,   0.1253,   0.12478,   0.12388,   0.12323,   0.12205,   0.12181,   0.12159,   0.12144,   0.12085,   0.12098,   0.12048,   0.12046,   0.12019,   0.11992,   0.11863,   0.11848,   0.11745,   0.11711,   0.1144,   0.11353,   0.1099,   0.11062,   0.10781,   0.10679,   0.10639,   0.10408};
  
  /*
  //make 425 nm to be mult mm
  double abs_length[nentries];
  for (unsigned i = 0; i < nentries; ++i)
  {
    abs_length[i] = Normabs_length[i] * mult * mm / 12.09;
  }
  */
  static const unsigned nentries = 2;
  static double photonE[nentries] = {1* eV, 6*eV};
  static double Normabs_length[nentries] = {1* mm, 1* mm};
  //flat abs
  double abs_length[nentries];
  for (unsigned i = 0; i < nentries; ++i)
  {
    abs_length[i] = Normabs_length[i] * mult;
  }

  G4MaterialPropertiesTable *tabler = material->GetMaterialPropertiesTable();
  tabler->RemoveProperty("ABSLENGTH");
  tabler->AddProperty("ABSLENGTH", photonE, abs_length, nentries);
}


//wls
G4Material* Make_Y11(){
 
  //--------------------------------------------------
  // WLSfiber Y11
  //--------------------------------------------------

  G4double density;
  std::vector<G4int> natoms;
  std::vector<G4String> elements;
  elements.push_back("C");
  natoms.push_back(5);
  elements.push_back("H");
  natoms.push_back(8);
  elements.push_back("O");
  natoms.push_back(2);

  density = 1.190 * g / cm3;

  G4Material* fPMMA = G4NistManager::Instance()->ConstructNewMaterial("PMMA", elements, natoms, density);

/*
  G4double photonEnergy[] =
  {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
   2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
   2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
   2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
   2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
   2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
   2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
   3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
   3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
   3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

  G4double absWLSfiber[] =
  {5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
   5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
   5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,1.10*m,
   1.10*m,1.10*m,1.10*m,1.10*m,1.10*m,1.10*m, 1.*mm, 1.*mm, 1.*mm, 1.*mm,
    1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm};

  G4double emissionFib[] =
  {0.05, 0.10, 0.30, 0.50, 0.75, 1.00, 1.50, 1.85, 2.30, 2.75,
   3.25, 3.80, 4.50, 5.20, 6.00, 7.00, 8.50, 9.50, 11.1, 12.4,
   12.9, 13.0, 12.8, 12.3, 11.1, 11.0, 12.0, 11.0, 17.0, 16.9,
   15.0, 9.00, 2.50, 1.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00,
   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
*/
  //Y11 from
  //http://kuraraypsf.jp/psf/ws.html
  /*
  G4double photonEnergy[] ={
2.862439029	*eV,2.840778923	*eV,2.815225831	*eV,2.788754124	*eV,2.756033268	*eV,2.741312912	*eV,2.729395813	*eV,2.721518432	*eV,2.707216697	*eV,2.707181235	*eV,2.693076188	*eV,2.693035246	*eV,2.67778668	*eV,2.669025741	*eV,2.668951059	*eV,2.665291033	*eV,2.66035048	*eV,2.660327649	*eV,2.659175227	*eV,2.659141012	*eV,2.651703078	*eV,2.645581665	*eV,2.645536511	*eV,2.640663232	*eV,2.63824174	*eV,2.634587028	*eV,2.632171074	*eV,2.626078227	*eV,2.624888601	*eV,2.616425667	*eV,2.616398064	*eV,2.609196999	*eV,2.599591613	*eV,2.59122036	*eV,2.585261508	*eV,2.584054715	*eV,2.578123376	*eV,2.57219786	*eV,2.560459892	*eV,2.553442102	*eV,2.541874304	*eV,2.53611406	*eV,2.530385027	*eV,2.516790443	*eV,2.507816703	*eV,2.494493004	*eV,2.482417875	*eV,2.469386576	*eV,2.457557584	*eV,2.447961085	*eV,2.438420061	*eV,2.42893409	*eV,2.418492889	*eV,2.417422599	*eV,2.408094304	*eV,2.397807939	*eV,2.394710008	*eV,2.386565943	*eV,2.37846796	*eV,2.376439519	*eV,2.370420214	*eV,2.366412214	*eV,2.360430055	*eV,2.354482536	*eV,2.348564913	*eV,2.343655498	*eV,2.331984293	*eV,2.322358317	*eV,2.312811483	*eV,2.301457528	*eV,2.293013889	*eV,2.27724336	*eV,2.266230965	*eV,2.259854057	*eV,2.249910639	*eV,2.243625119	*eV,2.232938727	*eV,2.221477748	*eV,2.213613724	*eV,2.203206732	*eV,2.188647275	*eV,2.182698942	*eV,2.174263731	*eV,2.161739991	*eV,2.146903611	*eV,2.130661072	*eV,2.114658849	*eV,2.091108237	*eV,2.081064162	*eV,2.059000322	*eV,2.041098506	*eV,2.032993517	*eV,2.019869653	*eV,2.006917391	*eV,1.992012672	*eV,1.962177447	*eV,1.948600936	*eV};

  G4double emissionFib[] = {0.0063, 0.00248, 0.0025, 0.00638, 0.01987, 0.04296, 0.07182, 0.09491, 0.16223, 0.13146, 
  0.2334, 0.20455, 0.28919, 0.38535, 0.32573, 0.42766, 0.50267, 0.48151, 0.56613, 0.53728, 
  0.59691, 0.67768, 0.64499, 0.71423, 0.75654, 0.79308, 0.82962, 0.86424, 0.89309, 0.95079, 
  0.92579, 0.9835, 1.00658, 1.00467, 0.99314, 0.97391, 0.96045, 0.92585, 0.88355, 0.83741, 
  0.79511, 0.76243, 0.73551, 0.70668, 0.69708, 0.6971, 0.71442, 0.73751, 0.75676, 0.76639, 
  0.76448, 0.73949, 0.72604, 0.6972, 0.67029, 0.63184, 0.59338, 0.56263, 0.5261, 0.50687, 
  0.4838, 0.46073, 0.42805, 0.39536, 0.37037, 0.34346, 0.31463, 0.30118, 0.28389, 0.26659, 
  0.25315, 0.23971, 0.21665, 0.2032, 0.18783, 0.17245, 0.15516, 0.13787, 0.12635, 0.10713, 
  0.09562, 0.08024, 0.06872, 0.0572, 0.04761, 0.03995, 0.03228, 0.02848, 0.02658, 0.01892, 
  0.01704, 0.01513, 0.00938, 0.00749, 0.00559, 0.00181, 0.0};
*/


//from paper
  G4double photonEnergy[] ={2.614246802	*eV,2.601505936	*eV,2.591177042	*eV,2.582063819	*eV,2.569633972	*eV,2.559556123	*eV,2.549551772	*eV,2.540728492	*eV,2.535242944	*eV,2.528692503	*eV,2.522175824	*eV,2.517853408	*eV,2.516775118	*eV,2.510319641	*eV,2.508176859	*eV,2.504964476	*eV,2.503897195	*eV,2.500700801	*eV,2.497507528	*eV,2.495386553	*eV,2.494327417	*eV,2.492206829	*eV,2.487986436	*eV,2.486933571	*eV,2.485876612	*eV,2.484825531	*eV,2.48063008	*eV,2.47958343	*eV,2.478532709	*eV,2.476443827	*eV,2.475400705	*eV,2.473317097	*eV,2.469155469	*eV,2.466047105	*eV,2.461909886	*eV,2.459848918	*eV,2.457786525	*eV,2.456759059	*eV,2.453681809	*eV,2.451634588	*eV,2.446526623	*eV,2.443470122	*eV,2.441439898	*eV,2.440426051	*eV,2.436374282	*eV,2.43133441	*eV,2.425309276	*eV,2.420310307	*eV,2.41433962	*eV,2.410373938	*eV,2.405440952	*eV,2.397580764	*eV,2.388804981	*eV,2.377200348	*eV,2.362854404	*eV,2.35244569	*eV,2.344927487	*eV,2.339322502	*eV,2.33374864	*eV,2.330043087	*eV,2.324508947	*eV,2.319001034	*eV,2.317172302	*eV,2.313519161	*eV,2.311699062	*eV,2.308063144	*eV,2.305346918	*eV,2.301730939	*eV,2.297232206	*eV,2.295433383	*eV,2.290062238	*eV,2.286498228	*eV,2.282054585	*eV,2.277628181	*eV,2.271461649	*eV,2.266201973	*eV,2.2609666	*eV,2.252297252	*eV,2.245408682	*eV,2.236002892	*eV,2.224147425	*eV,2.210752458	*eV,2.199162545	*eV,2.190139427	*eV,2.182814211	*eV,2.169907236	*eV,2.157947861	*eV,2.146119591	*eV,2.135195775	*eV,2.122077651	*eV,2.1060963	*eV,2.08886785	*eV,2.076311132	*eV,2.066084719	*eV,2.053803014	*eV,2.040245488	*eV,2.031072128	*eV,2.017119488	*eV,2.005408133	*eV,1.997226432	*eV,};
  G4double emissionFib[] = {0	,0.002538406	,0.005076812	,0.007615218	,0.017764443	,0.022841255	,0.027918067	,0.035533285	,0.048225315	,0.07106657	,0.096446232	,0.121825893	,0.134517923	,0.164974396	,0.187815651	,0.215733718	,0.241117779	,0.291877101	,0.317256762	,0.33756401	,0.362943671	,0.375635701	,0.401015362	,0.434010242	,0.459389903	,0.477158745	,0.515226037	,0.543148503	,0.560912946	,0.60406145	,0.626902705	,0.637056329	,0.680200433	,0.741117779	,0.786800289	,0.814718356	,0.829948792	,0.847713235	,0.865482077	,0.888323332	,0.916241399	,0.939087054	,0.959389903	,0.972081933	,0.984769564	,1	,0.984769564	,0.959389903	,0.94162106	,0.926395023	,0.908626181	,0.888323332	,0.875631302	,0.860405265	,0.840102416	,0.832487198	,0.81217995	,0.796953913	,0.774112658	,0.753805409	,0.72081053	,0.695430869	,0.672589614	,0.64974396	,0.626902705	,0.596446232	,0.573604977	,0.550759322	,0.525379661	,0.497461594	,0.467005121	,0.444163866	,0.423856617	,0.39593855	,0.368020483	,0.345179228	,0.317256762	,0.294415507	,0.271574252	,0.248728597	,0.228425748	,0.203046087	,0.185277245	,0.17512802	,0.16243599	,0.14974396	,0.137056329	,0.121825893	,0.114215074	,0.098984638	,0.081220195	,0.068528164	,0.060912946	,0.048225315	,0.045686909	,0.027918067	,0.025379661	,0.017764443	,0.002538406	,0	,};
  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);
  //std::cout<<"G4int nEntries:"<<nEntries<<std::endl;
  //for(int i=0;i<nEntries;i++){
  //  absWLSfiber[i]/=100;
  //}

  G4double energySmall[]= { 1.0 * eV, 6.0 * eV };
  G4double refractiveIndexWLSfiber[] = { 1.764, 1.764 }; //to make c/n = 1.7
  //G4double refractiveIndexWLSfiber[] = { 1.875, 1.875 }; //to make c/n = 1.6
  //cut at 500nm
  //G4double photonEnergy2[] = {1*eV, 2.4799*eV, 2.48*eV, 2.4801*eV, 6.0 * eV };
  //cut at emission end
  G4double photonEnergy2[] = {1*eV, 2.62*eV, 2.621*eV, 6.0 * eV };
  G4double absWLSfiber[] = {500.0*m, 500.0*m, 0.1*mm, 0.1*mm};
  G4double realabs_e[]={2.469788905	*eV,2.46085797	*eV,2.450328967	*eV,2.441473547	*eV,2.431120857	*eV,2.422141814	*eV,2.4119624	*eV,2.403099844	*eV,2.393371904	*eV,2.384929095	*eV,2.375382911	*eV,2.365587418	*eV,2.356763874	*eV,2.348034364	*eV,2.339321619	*eV,2.321568287	*eV,2.295443157	*eV,2.278841635	*eV,2.253935662	*eV,2.2134774	*eV,2.174693096	*eV,};
  G4double realabsWLSfiber[] = { 2.64495	*m,2.40574	*m,2.80068	*m,3.07617	*m,3.37919	*m,3.26865	*m,3.48896	*m,3.58062	*m,3.87446	*m,4.09481	*m,4.09454	*m,4.3516	*m,4.48921	*m,4.37866	*m,4.68171	*m,4.8742	*m,5.11237	*m,5.78281	*m,5.98422	*m,5.97371	*m,5.98159	*m,}; 
  const G4int nEntries_realabs = sizeof(realabs_e)/sizeof(G4double);
  //for(int i=0;i<nEntries_realabs;i++) realabs_e[i]=realabs_e[i];

  // Add entries into properties table
  G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();
  mptWLSfiber->AddProperty("RINDEX", energySmall, refractiveIndexWLSfiber, 2);
  mptWLSfiber->AddProperty("WLSCOMPONENT", photonEnergy, emissionFib, nEntries); //WLSCOMPONENT is the relative emission spectrum of the material as a function of the photon's momentum
  //mptWLSfiber->AddProperty("WLSABSLENGTH", photonEnergy, absWLSfiber, nEntries); //WLSABSLENGTH is the absorption length of the material as a function of the photon's momentum. -> re-emit
  mptWLSfiber->AddProperty("WLSABSLENGTH", photonEnergy2, absWLSfiber, 4); //WLSABSLENGTH is the absorption length of the material as a function of the photon's momentum. -> re-emit
  //mptWLSfiber->AddProperty("ABSLENGTH", realabs_e, realabsWLSfiber,nEntries_realabs); //ABSLENGTH is the absorption length of the material as a function of the photon's momentum. -> real absorption
  mptWLSfiber->AddConstProperty("WLSTIMECONSTANT", 11.5 * ns); //accounts for any time delay which may occur between absorption and re-emission of the photon, defalt delta

  fPMMA->SetMaterialPropertiesTable(mptWLSfiber);

  elements.clear();
  natoms.clear();
  return fPMMA;


}

G4Material* Make_Pethylene(){

  G4double density;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;

  //--------------------------------------------------
  // Cladding (polyethylene)
  //--------------------------------------------------

  elements.push_back("C");
  natoms.push_back(2);
  elements.push_back("H");
  natoms.push_back(4);

  density = 1.200 * g / cm3;

  G4Material* fPethylene =
    G4NistManager::Instance()->ConstructNewMaterial("Pethylene", elements, natoms, density);

  G4double refractiveIndexClad1[] = { 1.49, 1.49 };
  G4double energySmall[] = { 1.0 * eV, 6.0 * eV };

  G4double absClad[] = { 20.0 * m, 20.0 * m };

  // Add entries into properties table
  G4MaterialPropertiesTable* mptClad1 = new G4MaterialPropertiesTable();
  mptClad1->AddProperty("RINDEX", energySmall, refractiveIndexClad1, 2);
  mptClad1->AddProperty("ABSLENGTH", energySmall, absClad, 2);

  fPethylene->SetMaterialPropertiesTable(mptClad1);

  elements.clear();
  natoms.clear();
return fPethylene;
}

G4Material* Make_FPethylene(){

  G4double density;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;

  //--------------------------------------------------
  // Double Cladding (fluorinated polyethylene)
  //--------------------------------------------------

  elements.push_back("C");
  natoms.push_back(2);
  elements.push_back("H");
  natoms.push_back(4);

  density = 1.400 * g / cm3;

  G4Material* fFPethylene =
    G4NistManager::Instance()->ConstructNewMaterial("FPethylene", elements, natoms, density);

  G4double refractiveIndexClad2[] = { 1.42, 1.42 };
  G4double absClad[] = { 20.0 * m, 20.0 * m };

  // Add entries into properties table 
   G4double energySmall[] = { 2.0 * eV, 3.47 * eV };
  G4MaterialPropertiesTable* mptClad2 = new G4MaterialPropertiesTable();
  mptClad2->AddProperty("RINDEX", energySmall, refractiveIndexClad2, 2);
  mptClad2->AddProperty("ABSLENGTH", energySmall, absClad, 2);

  fFPethylene->SetMaterialPropertiesTable(mptClad2);
  elements.clear();
  natoms.clear();
return fFPethylene;
}
G4Material* Make_Polystyrene(){

  G4double density;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;

  //--------------------------------------------------
  // Polystyrene
  //--------------------------------------------------

  elements.push_back("C");
  natoms.push_back(8);
  elements.push_back("H");
  natoms.push_back(8);

  density = 1.050 * g / cm3;

  G4Material* fPolystyrene =
    G4NistManager::Instance()->ConstructNewMaterial("Polystyrene", elements, natoms, density);

  G4double energy[] = {
    2.00 * eV, 2.03 * eV, 2.06 * eV, 2.09 * eV, 2.12 * eV, 2.15 * eV, 2.18 * eV,
    2.21 * eV, 2.24 * eV, 2.27 * eV, 2.30 * eV, 2.33 * eV, 2.36 * eV, 2.39 * eV,
    2.42 * eV, 2.45 * eV, 2.48 * eV, 2.51 * eV, 2.54 * eV, 2.57 * eV, 2.60 * eV,
    2.63 * eV, 2.66 * eV, 2.69 * eV, 2.72 * eV, 2.75 * eV, 2.78 * eV, 2.81 * eV,
    2.84 * eV, 2.87 * eV, 2.90 * eV, 2.93 * eV, 2.96 * eV, 2.99 * eV, 3.02 * eV,
    3.05 * eV, 3.08 * eV, 3.11 * eV, 3.14 * eV, 3.17 * eV, 3.20 * eV, 3.23 * eV,
    3.26 * eV, 3.29 * eV, 3.32 * eV, 3.35 * eV, 3.38 * eV, 3.41 * eV, 3.44 * eV,
    3.47 * eV
  };
  G4double energySmall[] = { 2.0 * eV, 3.47 * eV };
  G4double refractiveIndexPS[] = { 1.50, 1.50 };

  G4double absPS[] = { 2. * cm, 2. * cm };

  G4double scintilFast[] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };

  // Add entries into properties table
  G4MaterialPropertiesTable* mptPolystyrene = new G4MaterialPropertiesTable();
  mptPolystyrene->AddProperty("RINDEX", energySmall, refractiveIndexPS, 2);
  mptPolystyrene->AddProperty("ABSLENGTH", energySmall, absPS, 2);
  mptPolystyrene->AddProperty("FASTCOMPONENT", energy, scintilFast, 50);
  mptPolystyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
  mptPolystyrene->AddConstProperty("RESOLUTIONSCALE", 1.0);
  mptPolystyrene->AddConstProperty("FASTTIMECONSTANT", 10. * ns);

  fPolystyrene->SetMaterialPropertiesTable(mptPolystyrene);

  // Set the Birks Constant for the Polystyrene scintillator
  fPolystyrene->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  elements.clear();
  natoms.clear();
return fPolystyrene;

}
G4Material* Make_Silicone(){

  G4double density;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;

  //--------------------------------------------------
  // Silicone (Template for Optical Grease)
  //--------------------------------------------------

  elements.push_back("C");
  natoms.push_back(2);
  elements.push_back("H");
  natoms.push_back(6);

  density = 1.060 * g / cm3;

  G4Material* fSilicone =
    G4NistManager::Instance()->ConstructNewMaterial("Silicone", elements, natoms, density);

  G4double refractiveIndexSilicone[] = { 1.46, 1.46 };

  G4double absClad[] = { 20.0 * m, 20.0 * m };
  G4double energySmall[] = { 2.0 * eV, 3.47 * eV };
  // Add entries into properties table
  G4MaterialPropertiesTable* mptSilicone = new G4MaterialPropertiesTable();
  mptSilicone->AddProperty("RINDEX", energySmall, refractiveIndexSilicone, 2);
  mptSilicone->AddProperty("ABSLENGTH", energySmall, absClad, 2);

  fSilicone->SetMaterialPropertiesTable(mptSilicone);

  elements.clear();
  natoms.clear();
  return fSilicone;

}
G4Material* Make_TiO2(){

  G4double density;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;

  //--------------------------------------------------
  // TiO2
  //--------------------------------------------------

  elements.push_back("Ti");
  natoms.push_back(1);
  elements.push_back("O");
  natoms.push_back(2);

  density = 4.26 * g / cm3;

  G4Material* TiO2 =
    G4NistManager::Instance()->ConstructNewMaterial("TiO2", elements, natoms, density);

  elements.clear();
  natoms.clear(); 
   return TiO2;

}

G4Material* Make_Coating(){

  G4double density;
  G4int ncomponents;
  G4double fractionmass;
  //--------------------------------------------------
  // Scintillator Coating - 15% TiO2 and 85% polystyrene by weight.
  //--------------------------------------------------

  density = 1.52 * g / cm3;

  G4Material* fCoating = new G4Material("Coating", density, ncomponents = 2);

  fCoating->AddMaterial(Make_TiO2(), fractionmass = 15 * perCent);
  fCoating->AddMaterial(Make_Polystyrene(), fractionmass = 85 * perCent);

   return fCoating;

}
