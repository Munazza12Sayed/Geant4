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
// * for the full disclaimer and the limitation of liabiliy.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEA4NT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//

#include "GasGapEventAction.hh"
#include "GasGapEventActionMessenger.hh"
#include "GasGapHit.h"
#include "GasGapTrajectory.hh"
#include "GasGapRunAction.hh"
#include "TrGEMAnalysis.hh"
#include "G4TrajectoryPoint.hh"   

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

GasGapEventAction::GasGapEventAction()
  //GasGapEventAction::GasGapEventAction(F02RunAction* run,HistoManager* histo)
  :drawFlag(false)
   //  :drawFlag(false),fRunAct(run),fHistoManager(histo)
{
  //  fPrintModulo = 100;
  new GasGapEventActionMessenger(this);
}

GasGapEventAction::~GasGapEventAction()
{;}

void GasGapEventAction::BeginOfEventAction(const G4Event* evt)
{

  //  G4int evtNb = evt->GetEventID();
  // if (evtNb%fPrintModulo == 0)
  //   G4cout << "\n---> Begin of event: " << evtNb << G4endl;

  // initialisation per event
  // fEnergyAbs = fEnergyGap = 0.;
  // fTrackLAbs = fTrackLGap = 0.;   

  
  TrGEMAnalysis::GetInstance()->PrepareNewEvent(evt) ;
  
  
  if(drawFlag)
    {
      G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
      if(pVVisManager)
        {
          G4UImanager::GetUIpointer()->ApplyCommand("/vis~/draw/current");
        }
    }
}

void GasGapEventAction::EndOfEventAction(const G4Event* evt )
{

  //  fRunAct->fillPerEvent(fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap);
  //fill histograms
  // //
  // fHistoManager->FillHisto(1, fEnergyAbs);
  // fHistoManager->FillHisto(2, fEnergyGap);
  // fHistoManager->FillHisto(3, fTrackLAbs);
  // fHistoManager->FillHisto(4, fTrackLGap);

  // fHistoManager->FillNtuple(fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap);


  TrGEMAnalysis::GetInstance()->EndOfEvent(evt) ;
  G4cout<<" End of event   "<<G4endl;

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  GasGapCollID    = SDman->GetCollectionID(colNam="GasGapHitCollection");
  
  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  GasGapHitCollection* GasGapHC    = 0;
  if(HCE)
    {
      GasGapHC    = (GasGapHitCollection*)(HCE->GetHC(GasGapCollID));
    }
  
  if(GasGapHC)
    {
      int n_hit = GasGapHC->entries();
      G4cout << "     " << n_hit
	     << " hits are stored in GasGap GasGapHitsCollection." << G4endl;
      G4double totE = 0;
      for(int i=0;i<n_hit;i++)
        { totE += (*GasGapHC)[i]->GetEdep(); }
      G4cout << "     Total energy deposition in GasGap : "
	     << totE / GeV << " (GeV)" << G4endl;
    }
  


  // get number of stored trajectories
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  G4cout << G4endl;
  G4cout<<" Number of trajectories in the event   "<<n_trajectories<<G4endl;
  G4cout << " Trajectories in tracker  "<<
    " ----------------------------------------------------"
	 << G4endl;
  
  
  for(G4int i=0; i<n_trajectories; i++)
    {

      GasGapTrajectory* trj =
	(GasGapTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);
      trj->ShowTrajectory();
      
      if(trj->GetParticleName()=="mu+" || trj->GetParticleName()=="mu-"){
      
      G4int trkid = trj->GetTrackID();
      G4int trjpoint = trj->GetPointEntries();
      
      G4cout<<" Number of points in trajectory   "<<trjpoint<<G4endl;

      G4double trajen = trj->GetInitialMomentum().mag();
      
      G4cout<<" Track initial momentum   "<<trajen<<G4endl;
      
      
      for(G4int  i=0 ; i < trjpoint; i++){
	
	G4TrajectoryPoint* aTrajectoryPoint = 
	  (G4TrajectoryPoint*)(trj->GetPoint(i));
	
	G4cout<<"  Loop points   "<<G4endl;
	
	G4double trajposx = aTrajectoryPoint->GetPosition().getX();
	G4double trajposy = aTrajectoryPoint->GetPosition().getY();
	G4double trajposz = aTrajectoryPoint->GetPosition().getZ();
	
	TrGEMAnalysis::GetInstance()->AddTrajPos(i,trajposx,trajposy,trajposz);
      }

      }
      G4cout<<"  Loop trajectories   "<<G4endl;
    }


  G4cout<<" out of  Loop trajectories   "<<G4endl;



  TrGEMAnalysis::GetInstance()->AddTrajInf(n_trajectories);

  G4cout<<" Add traj info   "<<G4endl;

  


  
  if(drawFlag)
    {
      G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
      if(pVVisManager)
        {
          if(GasGapHC)    GasGapHC->DrawAllHits();
          G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
        }
    }
}



