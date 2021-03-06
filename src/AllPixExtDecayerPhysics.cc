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
// $Id: AllPixExtDecayerPhysics.cc 72244 2013-07-12 08:49:56Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/src/AllPixExtDecayerPhysics.cc
/// \brief Implementation of the AllPixExtDecayerPhysics class
///
/// \author I. Hrivnacova; IPN, Orsay

#include "AllPixExtDecayerPhysics.hh"
#include "G4Pythia8Decayer.hh"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4Decay.hh>
#include <G4DecayWithSpin.hh>
#include <G4LambdaCDecayWithSpin.hh>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixExtDecayerPhysics::AllPixExtDecayerPhysics(const G4String& name)
  : G4VPhysicsConstructor(name)
{
/// Standard constructor
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixExtDecayerPhysics::~AllPixExtDecayerPhysics()
{
/// Destructor
}

//
// protected methods
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixExtDecayerPhysics::ConstructParticle()
{
/// Nothing to be done here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixExtDecayerPhysics::ConstructProcess()
{
/// Loop over all particles instantiated and add external decayer
/// to all decay processes if External decayer is set

  // Create Geant4 external decayer

    G4cout<<"!!!! Construct external physics!\n";
  G4Pythia8Decayer* extDecayer = new G4Pythia8Decayer();
  extDecayer->SetVerboseLevel(1);
     // The extDecayer will be deleted in G4Decay destructor

  aParticleIterator->reset();
  while ((*aParticleIterator)())
  {    
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if(particle->GetPDGEncoding() == 3122) particle->SetDecayTable((G4DecayTable*)0);
    
    if ( verboseLevel > 1 ) {
      G4cout << "!!! Setting ext decayer for: "
             <<  aParticleIterator->value()->GetParticleName() 
             << G4endl;
    } 

    G4ProcessVector* processVector = pmanager->GetProcessList();
    for (G4int i=0; i<processVector->length(); i++) {

      G4Decay* decay = dynamic_cast<G4Decay*>((*processVector)[i]);
      if ( decay ) decay->SetExtDecayer(extDecayer);
      if (decay && particle->GetPDGEncoding()==4122) pmanager->RemoveProcess(decay);
    }
    if(particle->GetPDGEncoding()==4122) {
        //G4DecayWithSpin* decayWithSpin = new G4DecayWithSpin();
        G4LambdaCDecayWithSpin* decayWithSpin = new G4LambdaCDecayWithSpin();
        decayWithSpin->SetExtDecayer(extDecayer);
        pmanager->AddProcess(decayWithSpin);
        // set ordering for PostStepDoIt and AtRestDoIt
        pmanager->SetProcessOrdering(decayWithSpin, idxAtRest);
        pmanager->SetProcessOrdering(decayWithSpin, idxPostStep);
        pmanager->DumpInfo();
    }

  }

  if ( verboseLevel > 0 ) {
    G4cout << "External decayer physics constructed." << G4endl;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
