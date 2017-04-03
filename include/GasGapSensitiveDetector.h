// $Id: SensitiveDetector.hh 22 2009-12-22 12:36:46Z schaelic $
 
#ifndef GasGapSensitiveDetector_h
#define GasGapSensitiveDetector_h 1

/**
 * @file
 * @brief Defines SensitiveDetector class.
 */

#include "G4VSensitiveDetector.hh"
//class DetectorConstruction;
//class RunAction;
 
#include "GasGapHit.h"              // <<- the hit "format" we define
class G4HCofThisEvent;
class G4ElectronIonPair;


// <<- means "H(it) C(ollections) of This Event"
 
//class G4Step;
//class G4TouchableHistory;

/*!
 * \brief Defines sensitve part of detector geometry.
 *
 * Stores Hits with 
 *  * deposited energy
 *  * position
 * in <i>Hit Collections of This Event</i>
 *
 * /sa ProcessHit()
 */
class GasGapSensitiveDetector : public G4VSensitiveDetector
{
 public:
  /// Constructor
  GasGapSensitiveDetector(G4String SDname);
  /// Destructor
  ~GasGapSensitiveDetector();
 
 public:
  /// @name methods from base class G4VSensitiveDetector
  //@{
  /// Mandatory base class method : it must to be overloaded:
  G4bool ProcessHits(G4Step *step, G4TouchableHistory *ROhist);
 
  /// (optional) method of base class G4VSensitiveDetector
  void Initialize(G4HCofThisEvent* HCE); 
   /// (optional) method of base class G4VSensitiveDetector
  void EndOfEvent(G4HCofThisEvent* HCE);
  //@}
   
 
 private:

  G4int fCluster;
  G4int fCluster2;
  //G4ThreeVector f;
  std::vector<G4ThreeVector>* fpos;
  G4ElectronIonPair* fElIonPair;
  

  typedef std::map<G4int,GasGapHit*> hitMap_t; //< Helper mapping layer number with hit
  hitMap_t hitMap;
  GasGapHitCollection*      hitCollection;
};

#endif

