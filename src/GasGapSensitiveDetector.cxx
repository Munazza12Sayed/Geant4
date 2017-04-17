#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4HCofThisEvent.hh"
#include "G4HCtable.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4LossTableManager.hh"
#include "G4ElectronIonPair.hh"
#include "G4UserSteppingAction.hh"

#include "TrGEMAnalysis.hh"
#include "GasGapSensitiveDetector.h"
#include "GasGapHit.h"

GasGapSensitiveDetector::GasGapSensitiveDetector(G4String SDname)
  : G4VSensitiveDetector(SDname)//,
    //    pdgCode(0),
    //    kineticEnergy(0.),
    //    positionX(0.),
    //    positionY(0.),
    //    positionZ(0.),
    //    momentumX(0.),
    //    momentumY(0.),
    //    momentumZ(0.),
    // charge(0),
    // neutSensitive(false),
    // elegap(false),
    // chargegap(false),
    // elegap_B(false),
    // chargegap_B(false),
    // primaryene(0.),
    // zinteraction(0),
    // contaPrimary(0),
    // contaInteraction(0),
    // contaSec(0),
    // contaSec_B(0),
    // contaTrack(0),
    // contaGar(0)
    
{
  G4cout << "*************************************" << G4endl ;
  G4cout << "*** Creating SD with name: "<< SDname << " ***" << G4endl ;  
  G4cout << "*************************************" << G4endl ;
  
  fElIonPair = G4LossTableManager::Instance()->ElectronIonPair();
  
  // Declaration of the hit collection name
  G4String myCollectionName = "GasGapHitCollection";
  collectionName.insert(myCollectionName);
  
  // ttTrack.clear();
  // ttTrack_B.clear();
  // ttTrack_Gar.clear();
  // postTrack.clear();
   
  //  zinteraction=-5;
   
  //  for(G4int t=0;t<9;t++){
   
  //    driftDep[t]=0.;
  //    driftDepI[t]=0.;
   
  //    transferDep[t]=0.;
  //    transferDepI[t]=0.;
  
  //    transfer2Dep[t]=0.;
  //    transfer2DepI[t]=0.;

  //    inductionDep[t]=0.;
  //    inductionDepI[t]=0.;
 
  //    driftDep_B[t]=0.;
  //    driftDepI_B[t]=0.;
   
  //    transferDep_B[t]=0.;
  //    transferDepI_B[t]=0.;
  
  //    transfer2Dep_B[t]=0.;
  //    transfer2DepI_B[t]=0.;

  //    inductionDep_B[t]=0.;
  //    inductionDepI_B[t]=0.;
  
  // }
}


GasGapSensitiveDetector::~GasGapSensitiveDetector()
{}

G4bool GasGapSensitiveDetector::ProcessHits(G4Step *step, G4TouchableHistory *)
{
  //This method is called every time a G4Step is performed in the logical volume
  //to which this SD is attached: the GAS GAP.
  
  //To identify where the step is we use the touchable navigation,
  //Remember we need to use PreStepPoint!
  
  // contaSec=0;
  // contaSec_B=0;
  // contaTrack=0;
  // contaGar=0;

  // G4double COINCIDENCE[9] = {1e10,50,100,500,1e3,1e4,1e5,1e6,1e8}; //50 ns time window fot the coincidence
   
  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
  G4int copyNo = touchable->GetVolume(0)->GetCopyNo();
  G4int layerIndex = copyNo;
  G4String volName = touchable->GetVolume(0)->GetName();

  //We get now the energy deposited by this step
  G4double edep = step->GetTotalEnergyDeposit(); 
  G4Track* track = step->GetTrack();
  G4double edepI= edep-step->GetNonIonizingEnergyDeposit(); 
  G4double edepNoI  = step->GetNonIonizingEnergyDeposit(); 
  
  TrGEMAnalysis::GetInstance()->AddEDepI(edepI);

  G4cout<<" Entering the SteppingAction    "<<G4endl;
  
  G4double energyprestep= step->GetPreStepPoint()->GetKineticEnergy();
  G4double energypoststep= step->GetPostStepPoint()->GetKineticEnergy();
  G4double deltae = step->GetDeltaEnergy();
  G4double SLength = step->GetStepLength();

  G4double e1Prim = std::sqrt(track->GetMomentum().mag2() + (105*MeV)*(105*MeV)) - 105*MeV;

  G4cout<<"  ========================================= "<<G4endl;
  
  G4cout<<" ======================================== "<<G4endl;
  G4cout<<" ======================================== "<<G4endl;
  G4cout<<" ======================================== "<<G4endl;
  G4cout<<"  Volume name                            "<<volName<<G4endl;
  G4cout<<"  Layer index                            "<<layerIndex<<G4endl;
  G4cout<<"  Step Length                     [mm]   "<<SLength<<G4endl;
  G4cout<<"  TotalEnergyDeposit              [MeV]  "<<edep<<G4endl;
  G4cout<<"  Energy deposited Non-Ionization [MeV]  "<<edepNoI<<G4endl;
  G4cout<<"  Energy deposited Ionization     [MeV]  "<<edepI<<G4endl;
  G4cout<<"  Kinetic energy pre-step         [MeV]  "<<energyprestep<<G4endl;
  G4cout<<"  Kinetic energy post-step        [MeV]  "<<energypoststep<<G4endl;
  G4cout<<"  Delta Energy                    [MeV]  "<<deltae<<G4endl;
  G4cout<<"  Delta Energy - Total            [MeV]  "<<fabs(deltae) - edep<<G4endl;
  G4cout<<"  (Delta Energy - Total) - EdepI  [MeV]  "<<(fabs(deltae) - edep) - edepI<<G4endl;
  G4cout<<"  Momentum magnitude              [MeV]  "<<track->GetMomentum()<<G4endl;
  G4cout<<" Kinetic energy formula           [MeV]  "<<e1Prim<<G4endl;
  
  double x= track->GetPosition().getX();
  double y= track->GetPosition().getY();
  double z= track->GetPosition().getZ();

  G4cout<<"  Track position   x [mm]  "<<x<<G4endl;
  G4cout<<"  Track position   y [mm]  "<<y<<G4endl;
  G4cout<<"  Track position   z [mm]  "<<z<<G4endl;

  TrGEMAnalysis::GetInstance()->AddnSteps();

  G4int nstep = TrGEMAnalysis::GetInstance()->GetnStep();

  G4cout<<" This is the Step number   "<<nstep<<G4endl;

  //  if( fabs(step->GetTrack()->GetParticleDefinition()->GetPDGEncoding()) == 13){
  TrGEMAnalysis::GetInstance()->AddEtotSteps(edepI);
  TrGEMAnalysis::GetInstance()->AdddeltaEtotSteps(deltae);
  //  }

  
  TrGEMAnalysis::GetInstance()->AddTrajPos_perstep(nstep,x,y,z);
  TrGEMAnalysis::GetInstance()->AddEDepI_perstep(nstep,edepI);
  TrGEMAnalysis::GetInstance()->AddETot_perstep(nstep,edep);
  TrGEMAnalysis::GetInstance()->AddEDelta_perstep(nstep,deltae);
  TrGEMAnalysis::GetInstance()->AddEDepSec_perstep(nstep,fabs(deltae)-edep);
  TrGEMAnalysis::GetInstance()->AddEDepImSec_perstep(nstep,(fabs(deltae)-edep)-edepI);
  TrGEMAnalysis::GetInstance()->AddSLength_perstep(nstep,SLength);


  if( fabs(step->GetTrack()->GetParticleDefinition()->GetPDGEncoding()) == 13){
    G4double edepmu = step->GetTotalEnergyDeposit();
    G4double edepImu = edepmu - step->GetNonIonizingEnergyDeposit();
    TrGEMAnalysis::GetInstance()->AddEDepImu_perstep(nstep,edepImu);
  }
  
  // double px= track->GetMomentum().getX();
  // double py= track->GetMomentum().getY();
  // double pz= track->GetMomentum().getZ();
  
  G4int pdg = track->GetParticleDefinition()->GetPDGEncoding();
  TrGEMAnalysis::GetInstance()->AddPDGID(nstep,pdg);
  G4StepPoint* point = step->GetPostStepPoint();
  const G4VProcess* proc = point->GetProcessDefinedStep();
  const G4String procname = proc->GetProcessName();
  G4int trackIndex=track->GetTrackID();
  double genz= track->GetVertexPosition().getZ(); 
  TrGEMAnalysis::GetInstance()->AddgenZ(nstep,genz);

  G4cout<<"  PDG ID               "<<pdg<<G4endl;
  G4cout<<"  Track Indx           "<<trackIndex<<G4endl;
  G4cout<<"  genZ          [mm]   "<<genz<<G4endl;
  G4cout<<"  parentID             "<<step->GetTrack()->GetParentID()<<G4endl;
  TrGEMAnalysis::GetInstance()->AddParentID(nstep,step->GetTrack()->GetParentID());
  G4cout<<"  KineticEnergy [MeV]  "<<step->GetTrack()->GetKineticEnergy()<<G4endl;

  fCluster = fElIonPair->SampleNumberOfIonsAlongStep(step);
  fCluster2 = fElIonPair->MeanNumberOfIonsAlongStep(step);
  fpos = fElIonPair->SampleIonsAlongStep(step);

 //  void DumpMeanEnergyPerIonPair();

  G4cout<<"  Number of Ion clusters  "<<fCluster<<G4endl;

  if(fCluster>2){

  TrGEMAnalysis::GetInstance()->AddNclust_perstep(nstep,fCluster);
  TrGEMAnalysis::GetInstance()->AddNclust2_perstep(nstep,fCluster2);
  
  
  G4int nionclust = fpos->size();
  G4cout<<"  Size of ion clusters along step  "<<nionclust<<G4endl;
  TrGEMAnalysis::GetInstance()->AddnIonAlStep(nstep,nionclust);
  
  for(int ni=0;ni<nionclust;ni++){

    G4ThreeVector temp; 
    G4cout<<" somelines    "<<G4endl;

    temp = (*fpos)[ni];

    G4cout<<" position in X    "<<temp.getX()<<G4endl;
    G4cout<<" position in Y    "<<temp.getY()<<G4endl;
    G4cout<<" position in Z    "<<temp.getZ()<<G4endl;

    TrGEMAnalysis::GetInstance()->AddnIonPosXStep(ni, temp.getX());
    TrGEMAnalysis::GetInstance()->AddnIonPosYStep(ni, temp.getY());
    TrGEMAnalysis::GetInstance()->AddnIonPosZStep(ni, temp.getZ());
      
    
  }

  }
  
  //  fpos->clear();
  //  const G4SteppingManager* pSM = fpSteppingManager;
  
  //secondaries                                                                                                                           
  //                                                                                                                                      
  //  const G4TrackVector* secondary = fpSteppingManager->GetSecondary();                                                                     
  // for (size_t lp=0; lp<(*secondary).size(); lp++) {                                                                                       
  //   particle = (*secondary)[lp]->GetDefinition();                                                                                         
  //   G4String name   = particle->GetParticleName();                                                                                        
  //   G4String type   = particle->GetParticleType();                                                                                        
  //   G4double charge = particle->GetPDGCharge();                                                                                           
  //   G4double energy = (*secondary)[lp]->GetKineticEnergy();                                                                               
  //   fRunAction->ParticleCount(name,energy);                                                                                               
  //   //energy spectrum                                                                                                                     
  //   if (charge > 3.)  ih = 2;                                                                                                             
  //   else if (particle == G4Gamma::Gamma())       ih = 3;                                                                                  
  //   else if (particle == G4Neutron::Neutron())   ih = 4;                                                                                  
  //   else if (particle == G4Proton::Proton())     ih = 5;                                                                                  
  //   else if (particle == G4Deuteron::Deuteron()) ih = 6;                                                                                  
  //   else if (particle == G4Alpha::Alpha())       ih = 7;                                                                                  
  //   else if (type == "nucleus")                  ih = 8;                                                                                  
  //   else if (type == "meson")                    ih = 9;                                                                                  
  //   else if (type == "baryon")                   ih = 10;                                                                                 
  //   analysis->FillH1(ih,energy);                                                                                                          
  //   //energy-momentum balance                                                                                                             
  //   G4ThreeVector momentum = (*secondary)[lp]->GetMomentum();                                                                             
  //   Q        += energy;                                                                                                                   
  //   Pbalance += momentum;                                                                                                                 
  //   //particle flag                                                                                                                       
  //   fParticleFlag[particle]++;                                                                                                            
  // }                                                  
  
  // G4String genprocess;
  // if(track->GetCreatorProcess()!=0) {
  //   const G4VProcess * genproc=track->GetCreatorProcess();  
  //   genprocess = genproc->GetProcessName();
  // }
  // else  {genprocess="primary";}
  // const G4LogicalVolume *  genLogVolume=track->GetLogicalVolumeAtVertex();
  // G4String genvolume= genLogVolume->GetName();
  
  // if ((*step->GetSecondary()).size()>0 && trackIndex==1 && contaInteraction == 0){
  //   zinteraction=z; 
  //   contaInteraction=1;
  // } 
  
   
   
  //     if( trackIndex==1 && (volName == "FakeBottom" || volName == "FakeTop")&& contaPrimary==0 ) {  
  //   //  if( trackIndex==1 ) {  
  //     	primaryene=energy;
  // 	 contaPrimary=1;
  // 	 neutSensitive=true;
  //  }

  //     // we're in drift gap
  //  if(volName == "GasGap1") {
  //     for(G4int t=0;t<9;t++){
  //        if(t==0) driftDep[t] += edep ;
  //        if(t==0) driftDepI[t] += edepI ;
  //        if(t>0 && track->GetGlobalTime()< COINCIDENCE[t]){
  //           driftDep[t] += edep ;
  //           driftDepI[t] += edepI ;
  //        }
  //     }
  //  }

  //     // we're in transfer1gap
  //  if(volName == "GasGap2") {
  //     for(G4int t=0;t<9;t++){
  //       if(t==0) transferDep[t] += edep ;
  //       if(t==0) transferDepI[t] += edepI ;
  //       if(t>0 && track->GetGlobalTime()< COINCIDENCE[t]){
  //        transferDep[t] += edep ;
  //        transferDepI[t] += edepI ;
  // 	 }
  //      } 
  //   }
 
  //    // we're in transfer2gap
  //   if(volName == "GasGap3") {
  //     for(G4int t=0;t<9;t++){
  //       if(t==0) transfer2Dep[t] += edep ;
  //       if(t==0) transfer2DepI[t] += edepI ;
  //       if(t>0 && track->GetGlobalTime()< COINCIDENCE[t]){
  //        transfer2Dep[t] += edep ;
  //        transfer2DepI[t] += edepI ;
  // 	 }
  //      } 
  //   }
  //    // we're in inductiongap
  //    if(volName == "GasGap4") {
  //     for(G4int t=0;t<9;t++){
  //       if(t==0) inductionDep[t] += edep ;
  //       if(t==0) inductionDepI[t] += edepI ;
  //       if(t>0 && track->GetGlobalTime()< COINCIDENCE[t]){
  //       // we're in induction1 gap
  //        inductionDep[t] += edep ;
  //        inductionDepI[t] += edepI ;
  // 	 }
  //      } 
  //   }
    
    
  //   //all for the second chamber in the super chamber
  //      // we're in drift gap
  //  if(volName == "GasGap1_B") {
  //     for(G4int t=0;t<9;t++){
  //        if(t==0) driftDep_B[t] += edep ;
  //        if(t==0) driftDepI_B[t] += edepI ;
  //        if(t>0 && track->GetGlobalTime()< COINCIDENCE[t]){
  //           driftDep_B[t] += edep ;
  //           driftDepI_B[t] += edepI ;
  //        }
  //     }
  //  }

  //     // we're in transfer1gap
  //  if(volName == "GasGap2_B") {
  //     for(G4int t=0;t<9;t++){
  //       if(t==0) transferDep_B[t] += edep ;
  //       if(t==0) transferDepI_B[t] += edepI ;
  //       if(t>0 && track->GetGlobalTime()< COINCIDENCE[t]){
  //        transferDep_B[t] += edep ;
  //        transferDepI_B[t] += edepI ;
  // 	 }
  //      } 
  //   }
 
  //    // we're in transfer2gap
  //   if(volName == "GasGap3_B") {
  //     for(G4int t=0;t<9;t++){
  //       if(t==0) transfer2Dep_B[t] += edep ;
  //       if(t==0) transfer2DepI_B[t] += edepI ;
  //       if(t>0 && track->GetGlobalTime()< COINCIDENCE[t]){
  //        transfer2Dep_B[t] += edep ;
  //        transfer2DepI_B[t] += edepI ;
  // 	 }
  //      } 
  //   }
  //    // we're in inductiongap
  //    if(volName == "GasGap4_B") {
  //     for(G4int t=0;t<9;t++){
  //       if(t==0) inductionDep_B[t] += edep ;
  //       if(t==0) inductionDepI_B[t] += edepI ;
  //       if(t>0 && track->GetGlobalTime()< COINCIDENCE[t]){
  //       // we're in induction1 gap
  //        inductionDep_B[t] += edep ;
  //        inductionDepI_B[t] += edepI ;
  // 	 }
  //      } 
  //   }
   
  //   //end second chamber in the super chamber
 
  
    
  //  // SAVE TRACKS IN LAYER 1
  //   //_________________________________________________________________________
  
  
  //    if(volName == "GasGap1" || volName == "GasGap2") {
  //     // we're in one of the gaps
     
  //    for(G4int T=0;T< ttTrack.size();T++){
  //       if (ttTrack[T]==trackIndex){
  // 	   contaSec=9999;
  // 	   break;
  // 	}
  //    }
    
  //    if(contaSec!=9999)  {
  //     // G4cout<<"New track trackIndex "<<trackIndex<<G4endl;
  //      TrGEMAnalysis::GetInstance()->SaveGapTrack(pdg, genprocess, genvolume,  genz, volName, energy );
  //       contaSec=trackIndex;
  // 	ttTrack.push_back(trackIndex);
  // 	// G4cout<<trackIndex <<" trackIndex pdg "<<pdg<<" genprocess "<<genprocess<<G4endl;
  
  //     }  
  
  //      //se la particella e'¨ un elettrone 
  //     if(pdg==11){
  //         elegap=true;
  //     }
  //      if(track->GetParticleDefinition()->GetPDGCharge()!=0){
  //         chargegap=true;
  //     }
      
      
  //     //FOR GARFIELD
  //     if(volName == "GasGap1"){
     
  //       for(G4int T=0;T< ttTrack_Gar.size();T++){
  //         if (ttTrack_Gar[T]==trackIndex){
  // 	     contaGar=9999;
  // 	     break;
  // 	  }
  //       }
  //     if(contaGar!=9999)  {
  //        // G4cout<<"New track trackIndex "<<trackIndex<<G4endl;
  //        TrGEMAnalysis::GetInstance()->SaveGarfieldQuantities(pdg, energy, x, y,  z, px, py, pz );
  //         contaGar=trackIndex;
  // 	  ttTrack_Gar.push_back(trackIndex);
  // 	  // G4cout<<trackIndex <<" trackIndex pdg "<<pdg<<" genprocess "<<genprocess<<G4endl;
  
  //     }  
     
  //     }//END GARFIELD
      
  // }


  //   // SAVE TRACKS IN LAYER 2
  //  //_________________________________________________________________________
  //    if(volName == "GasGap1_B" || volName == "GasGap2_B") {
  //     // we're in one of the gaps
      
  //    for(G4int T=0;T< ttTrack_B.size();T++){
  //       if (ttTrack_B[T]==trackIndex){
  // 	   contaSec_B=9999;
  // 	   break;
  // 	}
  //    }
  //    if(contaSec_B!=9999)  {
  //      // G4cout<<"New track B trackIndex "<<trackIndex<<G4endl;
  //      TrGEMAnalysis::GetInstance()->SaveGapTrack(pdg, genprocess, genvolume,  genz, volName, energy );
  //       contaSec_B=trackIndex;
  // 	ttTrack_B.push_back(trackIndex);
        
  // 	//stampa!!
  // 	// G4cout<<trackIndex <<" trackIndexpdg "<<pdg<<" genprocess "<<genprocess<<G4endl;

  
  //   }  
  //     //se la particella e'¨ un elettrone 
  //     if(pdg==11){
  //         elegap_B=true;
  //     }
  //      if(track->GetParticleDefinition()->GetPDGCharge()!=0){
  //         chargegap_B=true;
  //     }
  // }
  //   //_________________________________________________________________________


  //  if(volName == "FakeTop") {

  //       for(G4int T=0;T< postTrack.size();T++){
  //       if (postTrack[T]==trackIndex){
  // 	   contaTrack=9999;
  // 	   break;
  // 	}
  //    }
  //    if(contaTrack!=9999)  {
  //      // G4cout<<"POST track trackIndex "<<trackIndex<<G4endl;
  //       TrGEMAnalysis::GetInstance()->SavePostShieldTrack(pdg, energy );
  //       contaTrack=trackIndex;
  // 	postTrack.push_back(trackIndex);
  // 	//for(G4int aa=0;aa<postTrack.size();aa++){
  // 	//G4cout<<" "<<postTrack[aa];
  // 	//}
  // 	//G4cout<<G4endl;
  //   }  
  //  }
      
  //This line is used to store in Analysis class the energy deposited in this layer
  //The Analysis class will sum up this edep to the current event total energy in this layer
  //Pass the number directly to the Analysis manager. No Hits objects are created in 
  //this case
  //Analysis::GetInstance()->AddEDepHad(layerIndex,edep);

 
    
  
  
  //check if edep is from primary or secondary:
  // G4String isPri = step->GetTrack()->GetTrackID() == 1 ? "Yes" : "No";
  
  // // Tricks to implement hits 
  // hitMap_t::iterator it = hitMap.find(layerIndex) ;
  // GasGapHit* aHit = 0 ;
  // if( it != hitMap.end() ) {
  //   aHit = it->second ;
  // }
  // else {
  //   aHit = new GasGapHit(layerIndex) ;
  //   hitMap.insert( std::make_pair(layerIndex, aHit) ) ;
  //   hitCollection->insert(aHit) ;
  // }
  // aHit->AddEdep(edep) ;
  
  return true;
}

void GasGapSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
  hitCollection = new GasGapHitCollection(GetName(), collectionName[0]) ;
  static G4int HCID = -1 ;
  if(HCID<0) HCID = GetCollectionID(0) ; // <<-- this is to get an ID for collectionName[0]
  HCE->AddHitsCollection(HCID, hitCollection) ;

  G4cout<<"  +++++++++++++++++++++++++++++++++++++  "<<G4endl;
  
  G4cout<<"  Hit collection size                "<<hitCollection->entries()<<G4endl;
  
  // Reset map of hits
  hitMap.clear() ;
  //  fpos->clear();
}

void GasGapSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
  //  G4double ionizationPotential = 0.45*26*eV + 0.15*33*eV + 0.4*54*eV ; // Ar:CO2:CF4 (45:15:40)
  // Updated peer-reviewed values. Effective energy to generate a pair. (Sauli '77, Sharma)  
  // G4double ionizationPotential = 0.45*15.8*eV + 0.15*13.78*eV + 0.4*15.9*eV ; // Ar:CO2:CF4 (45:15:40)
  // These are values previously used. They represent the minimum ionization potential.
  //  G4int factor = 5 ;
  // TrGEMAnalysis::GetInstance()->SetSensitivity( driftDep,driftDepI, 									      transferDep,transferDepI,
  // 						transfer2Dep,transfer2DepI,
  // 						inductionDep,inductionDepI,
  // 						driftDep_B,driftDepI_B, 									      transferDep_B,transferDepI_B,
  // 						transfer2Dep_B,transfer2DepI_B,
  // 						inductionDep_B,inductionDepI_B) ;
  
  // TrGEMAnalysis::GetInstance()->SavePrimary(primaryene,zinteraction) ;
  
  // ttTrack.clear();
  // ttTrack_B.clear();
  // ttTrack_Gar.clear();
  // postTrack.clear();
       
  //     for(G4int t=0;t<9;t++){
        
  // 	driftDep[t] = 0. ;
  //       driftDepI[t] = 0. ;
     
  //       transferDep[t] = 0. ;
  //       transferDepI[t] = 0. ;
 
  //       transfer2Dep[t] = 0. ;
  //       transfer2DepI[t] = 0. ;

  //       inductionDep[t] = 0. ;
  //       inductionDepI[t] = 0. ;
 
  // 	driftDep_B[t] = 0. ;
  //       driftDepI_B[t] = 0. ;
     
  //       transferDep_B[t] = 0. ;
  //       transferDepI_B[t] = 0. ;
 
  //       transfer2Dep_B[t] = 0. ;
  //       transfer2DepI_B[t] = 0. ;

  //       inductionDep_B[t] = 0. ;
  //       inductionDepI_B[t] = 0. ;
  //    }

  //     //hitCollection->PrintAllHits() ;

  //     TrGEMAnalysis::GetInstance()->SetNeutronSensitivity(neutSensitive) ;
  //     TrGEMAnalysis::GetInstance()->SetEleGap(elegap, elegap_B) ;
  //     TrGEMAnalysis::GetInstance()->SetChargeGap(chargegap, chargegap_B) ;
  //     neutSensitive = false ;
  //     contaPrimary=0;
  //     contaInteraction=0;
  //     contaSec=0;
  //     contaSec_B=0;
  //     contaTrack=0;
  //     contaGar=0;
  //     elegap=false;
  //     chargegap=false;
  //     elegap_B=false;
  //     chargegap_B=false;
  //     primaryene=0.;
  //     zinteraction=-5.;
 
   G4cout<<"  +++++++++++++++++++++++++++++++++++++  "<<G4endl;

}

