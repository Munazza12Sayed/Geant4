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
#include "PrimaryGenerator.h"
#include "PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4HEPEvtInterface.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "TrGEMAnalysis.hh"

PrimaryGenerator::PrimaryGenerator()
{
  // const char* filename = "pythia_event.data";
  // HEPEvt = new G4HEPEvtInterface(filename);
  //create a messenger for this class
  
   time_t seed = time( NULL );
   fRandomEngine = new CLHEP::HepJamesRandom( static_cast< long >( seed ) );
   fRandomGauss = new CLHEP::RandGauss( fRandomEngine );
  
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
  SetDefaultKinematic();
  SetEnBeam(Enval);
  fSigmaPosition = 10.* mm;
  //  Enval=150*GeV;
  
  
  //create a messenger for this class  
  gunMessenger = new PrimaryGeneratorMessenger(this);

  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4String particleName;
  // G4ParticleDefinition* particle
  //   = particleTable->FindParticle(particleName="mu-");
  // fParticleGun->SetParticleDefinition(particle);
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  // fParticleGun->SetParticleEnergy(100.*keV);
  // fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,-2.4*cm));
  // particleGun = fParticleGun;
  
  //  messenger = new ExN04PrimaryGeneratorMessenger(this);
  //  useHEPEvt = false;
}


PrimaryGenerator::~PrimaryGenerator()
{
  //  delete HEPEvt;
  delete particleGun;
  delete gunMessenger;
}


void PrimaryGenerator::SetDefaultKinematic()
{
  //  G4int n_particle = 1;
  //  G4ParticleGun* fParticleGun = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
  = particleTable->FindParticle(particleName="mu-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.707,0.707));
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,-3.4*cm));
  particleGun = fParticleGun;

}
void PrimaryGenerator::SetEnBeam(G4double Enval){
  fParticleGun->SetParticleEnergy(Enval);
  G4cout<<" Energy of particle set to   "<<Enval<<"  [MeV]  "<<G4endl;
}

// void PrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
// {
//  if (fRndmBeam > 0.) {
//     G4double s_proj = fRndmBeam/std::sqrt(2.);
//     G4double x0, y0;     
//     do {
//         x0 = G4RandGauss::shoot(0.,s_proj);
//         y0 = G4RandGauss::shoot(0.,s_proj);      
//     } while ((x0*x0 + y0*y0) > fR2World);

// fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,fZ0));
  // if(useHEPEvt)
  // { HEPEvt->GeneratePrimaryVertex(anEvent); }
  // else
  // { particleGun->GeneratePrimaryVertex(anEvent); }

void PrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
{
      G4ThreeVector position = fParticleGun->GetParticlePosition();
    
     fRandomGauss->fire(0,20);
     fRandomGauss->fire(0,10);
     fParticleGun->SetParticlePosition(position);
     fParticleGun->GeneratePrimaryVertex(anEvent);
     
      //position.setY(dy);

  particleGun->GeneratePrimaryVertex(anEvent);
  
  TrGEMAnalysis::GetInstance()->AddPrimPos( fRandomGauss->fire(0,20),fRandomGauss->fire(0,10), particleGun->GetParticlePosition().getZ());

   G4cout<< "G4uniformra" <<G4UniformRand() << G4endl;
  G4cout<< "sigmapos" <<fSigmaPosition<< G4endl;

  //  TrGEMAnalysis::GetInstance()->AddPrimEn(particleGun->GetParticleGun()->GetParticleEnergy());

   G4cout<< "frand" <<  fRandomGauss->fire(0,20 )<< G4endl;
   G4cout<< "frand" <<  fRandomGauss->fire(0,10 )<< G4endl;
}


// void PrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
// {
//   // if(useHEPEvt)
//   // { HEPEvt->GeneratePrimaryVertex(anEvent); }
//   // else
//   // { particleGun->GeneratePrimaryVertex(anEvent); }
//   particleGun->GeneratePrimaryVertex(anEvent);
  
//   TrGEMAnalysis::GetInstance()->AddPrimPos(particleGun->GetParticlePosition().getX(),particleGun->GetParticlePosition().getY(), particleGun->GetParticlePosition().getZ());
//   //  TrGEMAnalysis::GetInstance()->AddPrimEn(particleGun->GetParticleGun()->GetParticleEnergy());
    
// }
