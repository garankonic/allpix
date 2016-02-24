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
    fMessenger(this),
    fVerboseLevel(0),
    fDecayProductsArray(0)
{
/// Standard constructor

  fDecayProductsArray = new ParticleVector();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Pythia8Decayer::~G4Pythia8Decayer()
{
/// Destructor

  delete fDecayProductsArray;
}

//
// private methods
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* G4Pythia8Decayer::
GetParticleDefinition(const Pythia8Particle* particle, G4bool warn) const
{
/// Return G4 particle definition for given TParticle

  // get particle definition from G4ParticleTable
  G4int pdgEncoding = particle->fKF;
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
G4Pythia8Decayer::CreateDynamicParticle(const Pythia8Particle* particle) const
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
                                   const Pythia8Particle* particle) const
{
/// Return particle vertex position.

  G4ThreeVector position 
     = G4ThreeVector(particle->fVx * cm,
                     particle->fVy * cm,
                     particle->fVz * cm);
  return position;
}                       
                        
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector G4Pythia8Decayer::GetParticleMomentum(
                                   const Pythia8Particle* particle) const
{
/// Return particle momentum.

  G4ThreeVector momentum 
     = G4ThreeVector(particle->fPx * GeV,
                     particle->fPy * GeV,
                     particle->fPz * GeV);
  return momentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4Pythia8Decayer::CountProducts(G4int channel, G4int particle)
{
/// Count number of decay products

   G4int np = 0;
   for ( G4int i=1; i<=5; i++ ) 
      if ( std::abs(Pythia8::Instance()->GetKFDP(channel,i) ) == particle )
        np++;
   return np;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
G4Pythia8Decayer::ForceParticleDecay(G4int particle, G4int product, G4int mult)
{
/// Force decay of particle into products with multiplicity mult

   Pythia8* Pythia8 = Pythia8::Instance();

   G4int kc =  Pythia8->Pycomp(particle);
   Pythia8->SetMDCY(kc,1,1);

   G4int ifirst = Pythia8->GetMDCY(kc,2);
   G4int ilast  = ifirst + Pythia8->GetMDCY(kc,3)-1;

   //
   //  Loop over decay channels
   for (G4int channel= ifirst; channel <= ilast; channel++) {
      if (CountProducts(channel,product) >= mult) {
         Pythia8->SetMDME(channel,1,1);
      } else {
         Pythia8->SetMDME(channel,1,0);
      }
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Pythia8Decayer::Decay(G4int pdg, const CLHEP::HepLorentzVector& p)
{
/// Decay a particle of type IDPART (PDG code) and momentum P.

   Pythia8::Instance()->Py1ent(0, pdg, p.e(), p.theta(), p.phi());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4Pythia8Decayer::ImportParticles(ParticleVector* particles)
{
/// Get the decay products into the passed PARTICLES vector

   return Pythia8::Instance()->ImportParticles(particles,"All");
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
  G4int nofParticles = ImportParticles(fDecayProductsArray);
  
  if ( fVerboseLevel > 0 ) {
    G4cout << "nofParticles: " <<  nofParticles << G4endl;
  }  

  // convert decay products Pythia8Particle type
  // to G4DecayProducts  
  G4DecayProducts* decayProducts
    = new G4DecayProducts(*(track.GetDynamicParticle()));

  G4int counter = 0;
  for (G4int i=0; i<nofParticles; i++) {

    // get particle from ParticleVector
    Pythia8Particle* particle = (*fDecayProductsArray)[i];
      
    G4int status = particle->fKS;
    G4int pdg = particle->fKF;
    if ( status>0 && status<11 && 
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
