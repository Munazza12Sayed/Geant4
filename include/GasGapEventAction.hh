//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//

#ifndef GasGapEventAction_h
#define GasGapEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "GasGapRunAction.hh"

class G4Event;
class GasGapRunAction;


class GasGapEventAction : public G4UserEventAction
{
  public:
  //  GasGapEventAction(F02RunAction*,HistoManager*);
  GasGapEventAction();
    virtual ~GasGapEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  // void AddAbs(G4double de, G4double dl) {fEnergyAbs += de; fTrackLAbs += dl;};
  // void AddGap(G4double de, G4double dl) {fEnergyGap += de; fTrackLGap += dl;};  
  
  private:

  //  F02RunAction*    fRunAct;
  //  HistoManager* fHistoManager;
  
  //  G4double  fEnergyAbs, fEnergyGap;
  // G4double  fTrackLAbs, fTrackLGap;
  
  // G4int     fPrintModulo;  
  
  G4int GasGapCollID;
  G4bool drawFlag;
  
public:
  inline void SetDrawFlag(G4bool val)
  { drawFlag = val; };
};

#endif

    
