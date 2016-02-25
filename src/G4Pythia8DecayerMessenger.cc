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
#include <G4UIcmdWithoutParameter.hh>

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Pythia8DecayerMessenger::G4Pythia8DecayerMessenger(
                               G4Pythia8Decayer* Pythia8Decayer)
  : G4UImessenger(),
    fPythia8Decayer(Pythia8Decayer),
    fDirectory(0),
    fVerboseCmd(0),
    fSetSeedCmd(0),
    fPythiaReadCmd(0),
    fPythiaReadFileCmd(0),
    fPythiaInitCmd(0)
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

  fPythiaReadCmd = new G4UIcommand("/Pythia8Decayer/Read",this);
  fPythiaReadCmd->SetGuidance("Read String");
  G4UIparameter* string = new G4UIparameter ("String to Read", 's', false);
  fPythiaReadCmd->SetParameter(string);

  fPythiaReadFileCmd = new G4UIcommand("/Pythia8Decayer/ReadFile",this);
  fPythiaReadFileCmd->SetGuidance("Read File");
  G4UIparameter* filename = new G4UIparameter ("File to Read", 's', false);
  fPythiaReadFileCmd->SetParameter(filename);

  fSetSeedCmd = new G4UIcmdWithAnInteger("/Pythia8Decayer/SetSeed",this);
  fSetSeedCmd->SetGuidance("Set initial seed");

  fPythiaInitCmd = new G4UIcmdWithoutParameter("/Pythia8Decayer/Init",this);
  fPythiaInitCmd->SetGuidance("Initialize Pythia with parameters");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Pythia8DecayerMessenger::~G4Pythia8DecayerMessenger()
{
/// Destructor

    delete fDirectory;
    delete fVerboseCmd;
    delete fPythiaReadCmd;
    delete fPythiaReadFileCmd;
    delete fPythiaInitCmd;
    delete fSetSeedCmd;
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
                ->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue)); \
    }
    else if(command == fSetSeedCmd) {
        G4int iseed = fSetSeedCmd->GetNewIntValue(newValue);
        fPythia8Decayer->GetPythia8()->SetSeed(iseed);
    }
    else if(command == fPythiaInitCmd) {
        fPythia8Decayer->Init();
    }
    else if(command == fPythiaReadCmd) {
        G4String s = newValue;
        fPythia8Decayer->GetPythia8()->ReadString(s.c_str());
    }
    else if(command == fPythiaReadFileCmd) {
        G4String s = newValue;
        fPythia8Decayer->GetPythia8()->ReadConfigFile(s.c_str());
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
