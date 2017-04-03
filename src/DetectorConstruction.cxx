#include "DetectorConstruction.h"

//Material manager
#include "G4NistManager.hh"
#include "G4SDManager.hh"
//Basic Units.
#include "G4SystemOfUnits.hh"

//Basic Volume and Placemen
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "GasGapSensitiveDetector.h"
#include "DetectorMessenger.hh"

//Types of Volumes
#include "G4Box.hh"

//Include necessary scorers
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4UserLimits.hh"

DetectorConstruction::DetectorConstruction() :
  fGasMat(0),fCuMat(0),fKAPTONMat(0),fFR4Mat(0),
  G4VUserDetectorConstruction()
{
  
  fDetectorMessenger = new DetectorMessenger(this);
}


DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

void DetectorConstruction::DefineMaterials() {
  
  G4NistManager* manager = G4NistManager::Instance() ;

  G4Element* elC  = manager->FindOrBuildElement(6);
  G4Element* elF  = manager->FindOrBuildElement(9);
  G4Element* elSi = manager->FindOrBuildElement(14);
  G4Element* elO  = manager->FindOrBuildElement(8);
  G4Element* elH  = manager->FindOrBuildElement(1);
 
 
  
  G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu") ;
  G4Material *KAPTON = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
  //G4Material *Fe= G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
  G4Material *Al= G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  fCuMat = Cu;
  fKAPTONMat = KAPTON;
  //fFeMat = Fe;
  fAlMat = Al;

  G4int numel(0), natoms(0) ;
  G4double density(0.), temperature(0.), pressure(0.), fractionMass(0.)  ;
  G4String name, symbol ;
  
  G4Material* SiO2 =  new G4Material("quartz",density= 2.200*g/cm3, numel=2);
  SiO2->AddElement(elSi, natoms=1);
  SiO2->AddElement(elO , natoms=2);
  
  //from http://www.physi.uni-heidelberg.de/~adler/TRD/TRDunterlagen/RadiatonLength/tgc2.htm
  //Epoxy (for FR4 )
  density = 1.2*g/cm3;
  G4Material* Epoxy = new G4Material("Epoxy" , density, numel=2);
  Epoxy->AddElement(elH, natoms=2);
  Epoxy->AddElement(elC, natoms=2);
  
  //FR4 (Glass + Epoxy)
  density = 1.86*g/cm3; 
  G4Material* FR4 = new G4Material("FR4"  , density, numel=2);
  FR4->AddMaterial(Epoxy, fractionMass=0.472);
  FR4->AddMaterial(SiO2, fractionMass=0.528);
  fFR4Mat = FR4;
  
  // gases at STP conditions 
  G4Material* Argon = manager->FindOrBuildMaterial("G4_Ar");
  G4Material* CarbonDioxide = manager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4Material* empty = manager->FindOrBuildMaterial("G4_Galactic");
  
  // Ar:CO2 (70:30) @ STP conditions
  G4double mixtureDensity = (Argon->GetDensity() * (70/100.0) + CarbonDioxide->GetDensity() * (30/100.0)) ;
  G4Material *ArCO2 = new G4Material("Ar/CO2",mixtureDensity,2) ;
  ArCO2->AddMaterial(Argon, 0.7) ;
  ArCO2->AddMaterial(CarbonDioxide, 0.3);

  // G4int ncomponents;
  // G4double fractionmass;
  
  // G4Material* ArCO2 =
  //   new G4Material("Ar/CO2",   density= 1000.8223*mg/cm3, ncomponents=2);
  // ArCO2->AddMaterial(Argon,  fractionmass=0.7844);
  // ArCO2->AddMaterial(CarbonDioxide, fractionmass=0.2156);

  
  fGasMat = ArCO2;
  //fGasMat=Fe;
  //fGasMat=Al;
  //  fGasMat =Cu;
  
}


/**This routine is used to construct the geometry to be simulated. This includes the
 * necessary materials to produce objects.
 */
G4VPhysicalVolume *DetectorConstruction::Construct()
{
  DefineMaterials() ;

  // SD Manager 
  G4SDManager* sdman = G4SDManager::GetSDMpointer() ;
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = false;
  
  //--------------Register Materials-------------------------------//
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  
  // Construct the World
  G4double world_sizeXY = 10*cm;
  G4double world_sizeZ  = 5*cm;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
	      0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);    //its size
  
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
			air,           //its material
			"World");      //its name
      
  G4VPhysicalVolume* physWorld = new G4PVPlacement(
						   0,    		//no rotation
						   G4ThreeVector(),	//at (0,0,0)
						   logicWorld,	//contataing logical volume
						   "World",	        //name
						   0,		//no mother volume
						   false,		//no boolean operation
						   0,		//copy number
						   checkOverlaps);	//overlap checking
      
  //--------------Build the required geometry here-----------------//

  //	G4Material* MatGasGap;
  //	double* GasGapThick=3.0*mm;


  // G4Box* solidDriftBoard = new G4Box("DrftBoard",4.*cm,4.*cm,0.5*cm);
  // G4LogicalVolume* logicDriftBoard = new G4LogicalVolume (solidDriftBoard,fFR4Mat,"DriftBoard");
  // G4ThreeVector posDriftBoard = G4ThreeVector(0,0,-1.5*cm);
  // new G4PVPlacement(0,posDriftBoard,logicDriftBoard,"DriftBoard",logicWorld,false,0,checkOverlaps);


  G4double gasgapwidth=1*cm;
  
  G4Box* solidGasGap = new G4Box("GasGap",4.*cm,4.*cm,gasgapwidth/2.);
  G4LogicalVolume* logicGasGap = new G4LogicalVolume (solidGasGap,fGasMat,"GasGap");
  G4ThreeVector posGasGap = G4ThreeVector(0,0,0);
  new G4PVPlacement(0,posGasGap,logicGasGap,"GasGap",logicWorld,false,0,checkOverlaps);
      
  G4String GasGapSDname="GasGap";
  GasGapSensitiveDetector* GasGapSD = new GasGapSensitiveDetector(GasGapSDname) ;
  sdman->AddNewDetector(GasGapSD) ;
  logicGasGap->SetSensitiveDetector(GasGapSD);
  
  //  G4double maxStep = 0.1*gasgapwidth;
  //  stepLimit = new G4UserLimits(maxStep);
  //  logicGasGap->SetUserLimits(stepLimit);

  // G4double ionizationPotential = 0.7*26*eV + 0.3*33*eV = 28.1 ; // Ar:CO2(70:30)
  if(0.0 == fGasMat->GetIonisation()->GetMeanEnergyPerIonPair()) {
    SetPairEnergy(28.1*eV);
  }

  //Return the world
  return physWorld;
}

/**In this routine we construct the Senditive Detector (SD). This is where the scorers are registered.
 */
void DetectorConstruction::ConstructSDandField() 
{
}


void DetectorConstruction::SetPairEnergy(G4double val)
{
  if(val > 0.0) {
    fGasMat->GetIonisation()->SetMeanEnergyPerIonPair(val);
  }
}                      
