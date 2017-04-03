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
/// \file field/field02/src/GasGapRunAction.cc
/// \brief Implementation of the GasGapRunAction class
//
// $Id$
// 

#include "GasGapRunAction.hh"
#include "TrGEMAnalysis.hh"

//#include "GasGapRunMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
//#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//GasGapRunAction::GasGapRunAction(HistoManager* histo)
GasGapRunAction::GasGapRunAction()
 : G4UserRunAction()
   //   fMessenger(0),
   //   fSaveRndm(0)
{
  //  fMessenger = new GasGapRunMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GasGapRunAction::~GasGapRunAction()
{
  //  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GasGapRunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  TrGEMAnalysis::GetInstance()->PrepareNewRun(aRun) ;

  // save Rndm status
  //  if (fSaveRndm > 0)
  //  { 
    //    CLHEP::HepRandom::showEngineStatus();
    //    CLHEP::HepRandom::saveEngineStatus("beginOfRun.rndm");
  //  }  
  //  G4UImanager* uiManager = G4UImanager::GetUIpointer();
   
  //  G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();

  //  if (visManager) uiManager->ApplyCommand("/vis/scene/notifyHandlers");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GasGapRunAction::EndOfRunAction(const G4Run* aRun)
{

  
  TrGEMAnalysis::GetInstance()->EndOfRun(aRun) ;

  //  if (G4VVisManager::GetConcreteInstance())
  //  {
  //    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  //  }

  // save Rndm status

  // if (fSaveRndm == 1)
  // { 
  //   CLHEP::HepRandom::showEngineStatus();
  //   CLHEP::HepRandom::saveEngineStatus("endOfRun.rndm");
  // }     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
