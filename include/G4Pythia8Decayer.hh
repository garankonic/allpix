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
// $Id: G4Pythia8Decayer.hh 72244 2013-07-12 08:49:56Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/include/G4Pythia8Decayer.hh
/// \brief Definition of the G4Pythia8Decayer class
//
#ifndef G4_Pythia8_DECAYER_H
#define G4_Pythia8_DECAYER_H

#include "G4VExtDecayer.hh"
#include "G4Pythia8DecayerMessenger.hh"
#include "Pythia8/Pythia.h"
#include "G4Pythia8.hh"

#include "globals.hh"

class G4Track;
class G4DecayProducts;

/// Pythia8 decayer
///
/// Implements the G4VExtDecayer abstract class using the Pythia8 interface.
/// According to TPythia8Decayer class in Root:
/// http://root.cern.ch/
/// see http://root.cern.ch/root/License.html

class G4Pythia8Decayer : public G4VExtDecayer
{
  public:

    G4Pythia8Decayer();
    virtual ~G4Pythia8Decayer();

    virtual G4DecayProducts* ImportDecayProducts(const G4Track& track);
    
    void SetVerboseLevel(G4int verboseLevel) { fVerboseLevel =  verboseLevel; }
    void  Init();
    G4Pythia8* GetPythia8() { return fPythia8; }
    
  private:

    /// Not implemented
    G4Pythia8Decayer(const G4Pythia8Decayer& right);
    /// Not implemented
    G4Pythia8Decayer& operator=(const G4Pythia8Decayer& right);
    
    G4ParticleDefinition*
    GetParticleDefinition(const Pythia8::Particle* p,G4bool warn = true) const;
    G4DynamicParticle* CreateDynamicParticle(const Pythia8::Particle* p) const;
    G4ThreeVector GetParticlePosition(const Pythia8::Particle* particle) const;
    G4ThreeVector GetParticleMomentum(const Pythia8::Particle* particle) const;
    
    void  Decay(G4int pdg, const CLHEP::HepLorentzVector& p);
    void  AppendParticle(G4int pdg, const CLHEP::HepLorentzVector& p);
    void  ClearEvent();
    

    G4Pythia8DecayerMessenger fMessenger;  ///< command messenger
    G4int            fVerboseLevel;        ///< verbose level
    G4Pythia8*       fPythia8;
};

// ----------------------------------------------------------------------------

#endif
