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
// $Id: G4Pythia8DecayerMessenger.cc 72244 2013-07-12 08:49:56Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/src/G4Pythia8DecayerMessenger.cc
/// \brief Implementation of the G4Pythia8DecayerMessenger class

// ----------------------------------------------------------------------------
// Messenger class that defines commands for G4Pythia8Decayer.
//
// Implements command
// - /Pythia8Decayer/verbose [level]
// - /Pythia8Decayer/forceDecayType [decayType]
// ----------------------------------------------------------------------------

#include "G4Pythia8DecayerMessenger.hh"
#include "G4Pythia8Decayer.hh"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAnInteger.hh>

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Pythia8DecayerMessenger::G4Pythia8DecayerMessenger(
                               G4Pythia8Decayer* Pythia8Decayer)
  : G4UImessenger(),
    fPythia8Decayer(Pythia8Decayer),
    fDirectory(0),
    fVerboseCmd(0),
    fDecayTypeCmd(0)
{
/// Standard constructor

  fDirectory = new G4UIdirectory("/Pythia8Decayer/");
  fDirectory->SetGuidance("G4Pythia8Decayer control commands.");

  fVerboseCmd 
    = new G4UIcmdWithAnInteger("/Pythia8Decayer/verbose", this);
  fVerboseCmd->SetGuidance("Set Pythia8Decayer verbose level");
  fVerboseCmd->SetParameterName("VerboseLevel", false);
  fVerboseCmd->SetRange("VerboseLevel >= 0 && VerboseLevel <= 1");
  fVerboseCmd->AvailableForStates(G4State_Idle);

  fDecayTypeCmd 
    = new G4UIcmdWithAnInteger("/Pythia8Decayer/forceDecayType", this);
  fDecayTypeCmd->SetGuidance("Force the specified decay type");
  fDecayTypeCmd->SetParameterName("DecayType", false);
  /*std::ostringstream os;
  os << "DecayType >=  " << kSemiElectronic 
     << " && DecayType <= " << kMaxDecay;
  fDecayTypeCmd->SetRange(os.str().c_str());*/
  fDecayTypeCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Pythia8DecayerMessenger::~G4Pythia8DecayerMessenger()
{
/// Destructor

  delete fDirectory;
  delete fVerboseCmd;
  delete fDecayTypeCmd;
}

//
// public methods
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Pythia8DecayerMessenger::SetNewValue(G4UIcommand* command,
       G4String newValue)
{ 
/// Apply command to the associated object.

  if(command == fVerboseCmd) { 
    fPythia8Decayer
      ->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue)); 
  }   
  else if(command == fDecayTypeCmd) { 
    fPythia8Decayer
      ->ForceDecayType(EDecayType(fDecayTypeCmd->GetNewIntValue(newValue))); 
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
