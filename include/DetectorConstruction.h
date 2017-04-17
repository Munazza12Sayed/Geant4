/**
 *
 */

class G4LogicalVolume ;
class G4PhysicalVolume ;

#ifndef DETECTORCONSTRUCTION_H
#define DETECTORCONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4Trd.hh"

#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4FieldManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"


class  G4UserLimits;
class DetectorMessenger;


class DetectorConstruction : public G4VUserDetectorConstruction
{
 private:
  G4Material*        fGasMat;
  G4Material*        fCuMat;
  G4Material*        fKAPTONMat;
  G4Material*        fFR4Mat;
  G4Material*        fFeMat;
  G4Material*        fAlMat;
  G4Material*        fEmptyMat;
  G4Material*        fAirMat;
  G4double           tripleGemThinBase; 
  G4double           tripleGemLargeBase;
  G4double           tripleGemHeight;

  G4UserLimits* stepLimit;             // pointer to user step limits
  DetectorMessenger* fDetectorMessenger;
  //void DefineMaterials() ;
	       
  std::vector<G4Trd*>           trdCollection ;
  std::vector<G4LogicalVolume*> trdLogCollection ;

  G4ElectricField*        fEMfield;
  G4EqMagElectricField*   fEquation;
  G4MagIntegratorStepper* fStepper;
  G4FieldManager*         fFieldMgr;
  G4double                fMinStep ;
  G4ChordFinder*          fChordFinder ;
  G4MagInt_Driver*        fIntgrDriver;

		
 public:
  DetectorConstruction();
  virtual ~DetectorConstruction();
  void SetPairEnergy(G4double);		
		
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

  void DefineMaterials() ;
  //G4VPhysicalVolume* Construct() ;
  G4Trd* Trapezoid(G4String name, G4double width) ;
  void PlaceGeometry(G4RotationMatrix *pRot, G4ThreeVector tlate, G4LogicalVolume *pMotherLogical) ; 
	
};

#endif
