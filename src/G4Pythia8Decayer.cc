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
// $Id: G4Pythia8Decayer.cc 81443 2014-05-28 14:26:55Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/src/G4Pythia8Decayer.cc
/// \brief Implementation of the G4Pythia8Decayer class

// ----------------------------------------------------------------------------
// According to TPythia8Decayer class in Root:
// http://root.cern.ch/
// see http://root.cern.ch/root/License.html
// ----------------------------------------------------------------------------

#include "G4Pythia8Decayer.hh"
#include "Pythia8/Pythia.h"

#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include <CLHEP/Vector/LorentzVector.h>

#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Pythia8Decayer::G4Pythia8Decayer()
  : G4VExtDecayer("G4Pythia8Decayer"),
    fPythia8(new G4Pythia8()),
    fMessenger(this),
    fVerboseLevel(0)
{
/// Standard constructor

;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Pythia8Decayer::~G4Pythia8Decayer()
{
/// Destructor

  delete fPythia8;
}

//
// private methods
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* G4Pythia8Decayer::
GetParticleDefinition(const Pythia8::Particle* particle, G4bool warn) const
{
/// Return G4 particle definition for given TParticle

  // get particle definition from G4ParticleTable
  G4int pdgEncoding = particle->id();
  G4ParticleTable* particleTable 
    = G4ParticleTable::GetParticleTable();                
  G4ParticleDefinition* particleDefinition = 0;    
  if (pdgEncoding != 0) 
    particleDefinition = particleTable->FindParticle(pdgEncoding);

  if ( particleDefinition == 0 && warn) {
    G4cerr 
      << "G4Pythia8Decayer: GetParticleDefinition: " << std::endl
      << "G4ParticleTable::FindParticle() for particle with PDG = " 
      << pdgEncoding 
      << " failed." << std::endl;
  }
  
  return particleDefinition;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DynamicParticle*
G4Pythia8Decayer::CreateDynamicParticle(const Pythia8::Particle* particle) const
{ 
/// Create G4DynamicParticle.

  // get particle properties
  const G4ParticleDefinition* particleDefinition 
    = GetParticleDefinition(particle);    
  if ( ! particleDefinition ) return 0;  
        
  G4ThreeVector momentum = GetParticleMomentum(particle);

  // create G4DynamicParticle
  G4DynamicParticle* dynamicParticle 
    = new G4DynamicParticle(particleDefinition, momentum);
  
  return dynamicParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector G4Pythia8Decayer::GetParticlePosition(
                                   const Pythia8::Particle* particle) const
{
/// Return particle vertex position.

  G4ThreeVector position 
     = G4ThreeVector(particle->xProd() * mm,
                     particle->yProd() * mm,
                     particle->zProd() * mm);
  return position;
}                       
                        
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector G4Pythia8Decayer::GetParticleMomentum(
                                   const Pythia8::Particle* particle) const
{
/// Return particle momentum.

  G4ThreeVector momentum 
     = G4ThreeVector(particle->px() * GeV,
                     particle->py() * GeV,
                     particle->pz() * GeV);
  return momentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Pythia8Decayer::Init() {
    fPythia8->Pythia8()->init();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Pythia8Decayer::Decay(G4int pdg, const CLHEP::HepLorentzVector& p)
{
    ClearEvent();
    AppendParticle(pdg, p);
    G4int idPart = fPythia8->Pythia8()->event[0].id();
    //G4cout<<"Decaying: "<<idPart<<"\n";
    //fPythia8->Pythia8()->particleData.mayDecay(idPart,true);
    fPythia8->Pythia8()->moreDecays();
    if (fVerboseLevel > 0) fPythia8->EventListing();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Pythia8Decayer::AppendParticle(G4int pdg, const CLHEP::HepLorentzVector& p)
{
    fPythia8->Pythia8()->event.append(pdg, 11, 0, 0, p.px(), p.py(), p.pz(), p.e(), p.m());
}
void G4Pythia8Decayer::ClearEvent()
{
    fPythia8->Pythia8()->event.clear();
}

//
// public methods
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DecayProducts* G4Pythia8Decayer::ImportDecayProducts(const G4Track& track)
{
/// Import decay products

  // get particle momentum
  G4ThreeVector momentum = track.GetMomentum(); 
  G4double etot = track.GetDynamicParticle()->GetTotalEnergy();;  
  CLHEP::HepLorentzVector p;    
  p[0] = momentum.x() / GeV;
  p[1] = momentum.y() / GeV;
  p[2] = momentum.z() / GeV;
  p[3] = etot         / GeV;
  
  // get particle PDG
  // ask G4Pythia8Decayer to get PDG encoding
  // (in order to get PDG from extended TDatabasePDG
  // in case the standard PDG code is not defined)
  G4ParticleDefinition* particleDef = track.GetDefinition();
  G4int pdgEncoding = particleDef->GetPDGEncoding();

  // let Pythia8Decayer decay the particle
  // and import the decay products
  Decay(pdgEncoding, p);
  G4int nofParticles = fPythia8->GetN();
  
  if ( fVerboseLevel > 0 ) {
    G4cout << "nofParticles: " <<  nofParticles << G4endl;
  }  

  // convert decay products Pythia8::Particle type
  // to G4DecayProducts  
  G4DecayProducts* decayProducts
    = new G4DecayProducts(*(track.GetDynamicParticle()));

  G4int counter = 0;
  for (G4int i=1; i<=nofParticles; i++) {

    // get particle from ParticleVector
    Pythia8::Particle* particle = &(fPythia8->Pythia8()->event[i]);
      
    G4int status = particle->status();
    G4int mother = particle->mother1();
    G4int pdg = particle->id();
    if ( /*status>0 &&*/ mother==0 &&
         std::abs(pdg)!=12 && std::abs(pdg)!=14 && std::abs(pdg)!=16 ) {
      // pass to tracking final particles only;
      // skip neutrinos

      if ( fVerboseLevel > 0 ) {
        G4cout << "  " << i << "th particle PDG: " << pdg << "   ";
      }  
            
      // create G4DynamicParticle 
      G4DynamicParticle* dynamicParticle 
        = CreateDynamicParticle(particle);

      if (dynamicParticle) {

        if ( fVerboseLevel > 0 ) {
          G4cout << "  G4 particle name: " 
                 << dynamicParticle->GetDefinition()->GetParticleName()
                 << G4endl;
        }         

        // add dynamicParticle to decayProducts
        decayProducts->PushProducts(dynamicParticle);
        
        counter++;
      }
    }       
  }                             
  if ( fVerboseLevel > 0 ) {
    G4cout << "nofParticles for tracking: " <<  counter << G4endl;
  }  

  return decayProducts;
}
