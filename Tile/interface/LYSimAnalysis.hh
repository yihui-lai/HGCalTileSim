#ifndef LYSimAnalysis_HH_
#define LYSimAnalysis_HH_

#include "G4Event.hh"
#include "G4Run.hh"

#include <fstream>
#include <iostream>

#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"

class LYSimDetectorConstruction;
class LYSimPrimaryGeneratorAction;
class LYSimProtonGeneratorAction;
class LYSimScintillation;
class LYSimFormat;
class LYSimRunFormat;

class TFile;
class TTree;

class LYSimAnalysis
{
public:
  inline static LYSimAnalysis*
  GetInstance()
  {
    if( LYSimAnalysis::singleton == NULL ){
      LYSimAnalysis::singleton = new LYSimAnalysis();
    }
    return LYSimAnalysis::singleton;
  }

  virtual ~LYSimAnalysis();

  // Setting detector output
  void
  SetDetector( LYSimDetectorConstruction* det )
  { DetectorConstruction = det; };
  void
  SetGeneratorAction( LYSimPrimaryGeneratorAction* genaction )
  { generatorAction = genaction; };
  void SetProtonGeneratorAction( LYSimProtonGeneratorAction* genaction )
   { protonAction = genaction; }

  void
  SetOutputFile( const std::string& x ){ filename = x;}
  std::string
  GetOutputFile() const { return filename; }

  void PrepareExperiment();
  void PrepareNewRun( const G4Run* aRun );
  void PrepareNewEvent( const G4Event* anEvent );
  void EndOfEvent( const G4Event* anEvent );
  void EndOfRun( const G4Run* aRun );
  void EndOfExperiment();
  void addenergy(double, double);
  void addgenphoton();
  bool IsGenProton();
  void AddTracker(int, double, double, double);
  void addgenz(double, double);

  // Embedded class for RunAction handling
  class RunAction : public G4UserRunAction
  {
public:
    RunAction(){}
    ~RunAction(){}
    inline void
    BeginOfRunAction( const G4Run* run )
    { LYSimAnalysis::GetInstance()->PrepareNewRun( run ); }
    void
    EndOfRunAction( const G4Run* run )
    { LYSimAnalysis::GetInstance()->EndOfRun( run ); }
  };

  class EventAction : public G4UserEventAction
  {
public:
    EventAction(){}
    virtual ~EventAction(){}
    inline void
    BeginOfEventAction( const G4Event* event )
    { LYSimAnalysis::GetInstance()->PrepareNewEvent( event ); }
    inline void
    EndOfEventAction( const G4Event* event )
    { LYSimAnalysis::GetInstance()->EndOfEvent( event ); }
  };

private:
  // Singleton constructor.
  LYSimAnalysis();
  static LYSimAnalysis* singleton;

  LYSimDetectorConstruction* DetectorConstruction;
  LYSimPrimaryGeneratorAction* generatorAction;
  LYSimProtonGeneratorAction* protonAction;

public:
  std::string filename;
  TFile* file;
  TTree* tree;
  TTree* runtree;
  LYSimFormat* format;
  LYSimRunFormat* runformat;

  int GetNPhotons( const G4Event* event );
};

#endif/* LYSimAnalysis_HH_ */
