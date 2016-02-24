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
// $Id: G4Pythia8DecayerMessenger.hh 72244 2013-07-12 08:49:56Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/include/G4Pythia8DecayerMessenger.hh
/// \brief Definition of the G4Pythia8DecayerMessenger class

#ifndef G4_Pythia8_DECAYER_MESSENGER_H
#define G4_Pythia8_DECAYER_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class G4Pythia8Decayer;

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

/// Messenger class that defines commands for G4Pythia8Decayer.
///
/// Implements command
/// - /Pythia8Decayer/verbose [level]
/// - /Pythia8Decayer/forceDecayType [decayType]

class G4Pythia8DecayerMessenger : public G4UImessenger
{
  public:

    G4Pythia8DecayerMessenger(G4Pythia8Decayer* Pythia8Decayer);
    virtual ~G4Pythia8DecayerMessenger();
   
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  private:

    /// Not implemented
    G4Pythia8DecayerMessenger();
    /// Not implemented
    G4Pythia8DecayerMessenger(const G4Pythia8DecayerMessenger& right);
    /// Not implemented
    G4Pythia8DecayerMessenger& operator=(const G4Pythia8DecayerMessenger& r);

  private:

    G4Pythia8Decayer*      fPythia8Decayer;    ///< associated class
    G4UIdirectory*         fDirectory;         ///< command directory
    G4UIcmdWithAnInteger*  fVerboseCmd;        ///< command: verbose
    G4UIcmdWithAnInteger*  fDecayTypeCmd;      ///< command: forceDEcayeType
};

// ----------------------------------------------------------------------------

#endif
