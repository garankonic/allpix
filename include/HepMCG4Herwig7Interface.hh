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
/// \file eventgenerator/HepMC/HepMCEx03/include/HepMCG4Herwig7Interface.hh
/// \brief Definition of the HepMCG4Herwig7Interface class
//
//

#ifndef HEPMC_G4_Herwig7_INTERFACE_H
#define HEPMC_G4_Herwig7_INTERFACE_H

#include "HepMCG4Interface.hh"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"

/// A generic interface class with Herwig7 event generator via HepMC.

class HepMCG4Herwig7Interface : public HepMCG4Interface {
protected:
   G4int verbose;
  virtual HepMC::GenEvent* GenerateHepMCEvent();

private:
   typedef std::vector<std::string> CommandVector;
   CommandVector m_herwigCommandVector;
   /// ThePEG generator object.
   ThePEG::EGPtr m_hw;
   /// ThePEG event object.
   ThePEG::EventPtr m_event;
   double m_xsscale;

public:
  HepMCG4Herwig7Interface();
  ~HepMCG4Herwig7Interface();

  // set/get methods
  void SetVerboseLevel(G4int i);
  G4int GetVerboseLevel() const;

  // random numbers operations
  void PrintRandomStatus(std::ostream& ostr=G4cout) const;

  // setup user parameters (empty in default).
  // Implement your parameters in a delived class if you want.
  virtual void SetUserParameters();

  virtual void Print() const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void HepMCG4Herwig7Interface::SetVerboseLevel(G4int i)
{
  verbose= i;
}

inline G4int HepMCG4Herwig7Interface::GetVerboseLevel() const
{
  return verbose;
}

#endif
