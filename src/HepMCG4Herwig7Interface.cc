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
/// \file eventgenerator/HepMC/HepMCEx03/src/HepMCG4Herwig7Interface.cc
/// \brief Implementation of the HepMCG4PythiaInterface class for Herwig7
//

//#ifdef G4LIB_USE_Herwig7

#include "HepMCG4Herwig7Interface.hh"
#include "HepMC/GenEvent.h"

#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DynamicLoader.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/StandardSelectors.h"

#include <boost/algorithm/string.hpp>
#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

// Setup HepMC traits definition for ThePEG's converter to work
#include "ThePEG/Vectors/HepMCConverter.h"
namespace ThePEG {
  template<>
  struct HepMCTraits<HepMC::GenEvent>
    : public HepMCTraitsBase<HepMC::GenEvent,
                             HepMC::GenParticle,
                             HepMC::GenVertex,
                             HepMC::Polarization,
                             HepMC::PdfInfo>
  {
    static bool hasUnits() {
      #ifdef HEPMC_HAS_UNITS
      return true;
      #else
      return false;
      #endif
    }
  };
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4Herwig7Interface::HepMCG4Herwig7Interface()
   : verbose(0)
{
    ThePEG::Repository::exitOnError() = 1;
    ThePEG::Repository::load("/home/maren/GEANT4/Herwig-7.0.1-install/share/Herwig/defaults/HerwigDefaults.rpo");
    /// Prepend the defaults to the command vector, so that *everything* gets dumped
    CommandVector defaultcmds;
    defaultcmds.push_back("cd /Herwig/Generators");
    defaultcmds.push_back("set /Herwig/Generators/LHCGenerator:RandomNumberGenerator:Seed 12345");
    defaultcmds.push_back("set /Herwig/Generators/LHCGenerator:DebugLevel 1");
    defaultcmds.push_back("set /Herwig/Generators/LHCGenerator:PrintEvent 1");
    defaultcmds.push_back("set /Herwig/Generators/LHCGenerator:MaxErrors 1000000");
    defaultcmds.push_back("set /Herwig/Generators/LHCGenerator:NumberOfEvents 1000000000");
    defaultcmds.push_back("set /Herwig/Generators/LHCGenerator:UseStdout Yes");
    //defaultcmds.push_back("read LHC-Matchbox.in");

    CommandVector dirs;
    dirs.push_back("/home/maren/GEANT4/ThePEG-2.0.1-install/share");
    dirs.push_back("/home/maren/GEANT4/Herwig-7.0.1-install/share/Herwig");
    ThePEG::Repository::appendReadDir(dirs);

    foreach (const std::string& str, ThePEG::DynamicLoader::allPaths()) {
        G4cout<<str.c_str()<<G4endl;
    }
    defaultcmds.push_back("insert /Herwig/MatrixElements/SimpleQCD:MatrixElements[0] /Herwig/MatrixElements/MEHeavyQuark");
    defaultcmds.push_back("set /Herwig/MatrixElements/MEHeavyQuark:QuarkType Bottom");

    // Append config directives from job options to the defaults
    CommandVector cmds = defaultcmds;
    cmds.insert(cmds.begin()+cmds.size(), m_herwigCommandVector.begin(), m_herwigCommandVector.end());
    // Apply the config commands
    G4cout<<"Processing default and job option commands"<<G4endl;
    std::string commands = "";
    foreach (const std::string& cmd, cmds) {
      commands += "\n" + cmd;
      const size_t iNonSpace = cmd.find_first_not_of(" ");
      // Only run the command if it's not just whitespace or a comment
      if (iNonSpace != std::string::npos && cmd.data()[iNonSpace] != '#') {
        const std::string reply = ThePEG::Repository::exec(cmd, std::cout);
        if (!reply.empty()) {
          if (reply.find("Error") != std::string::npos) {
            G4cout<<"!!! Herwig++ error: "<<reply.c_str()<<G4endl;
          } else {
            G4cout<<"\tHerwig++ info: "<<reply.c_str()<<G4endl;
          }
        }
      }
    }
    G4cout<<"Updating repository"<<G4endl;
    ThePEG::Repository::update();

    // Dump out the config commands, with an extra saverun to make life easier
    std::ostringstream ss_cmds;
    ss_cmds << commands << "\n\n"
            << "# Extra saverun for standalone convenience: Athena doesn't execute this\n"
            << "saverun " <<  " /Herwig/Generators/LHCGenerator\n";
    const std::string dumpcommands = ss_cmds.str();
    if (1) {
      std::ofstream f("./dump.txt");
      f << dumpcommands;
      f.close();
    }

    // Make a "run" object from the config repository.
    ThePEG::EGPtr tmpEG = ThePEG::BaseRepository::GetObject<ThePEG::EGPtr>("/Herwig/Generators/LHCGenerator");
    try {
      m_hw = ThePEG::Repository::makeRun(tmpEG, "run000");
    } catch (ThePEG::Exception& e) {
      G4cout<<"!!! Exception in ThePEG: " << e.what();
      throw;
    } catch (std::exception& e) {
      G4cout<<"!!! STL exception: " << e.what();
      throw;
    }

    m_hw->initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4Herwig7Interface::~HepMCG4Herwig7Interface()
{
    G4cout << "MetaData: generator = Herwig++ " << G4endl;
    G4cout << std::scientific << std::setprecision(5) << "MetaData: cross-section (nb) = " << m_hw->eventHandler()->integratedXSec()*m_xsscale/ThePEG::nanobarn << G4endl;
    m_hw->finalize();
    ThePEG::Repository::cleanup();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Herwig7Interface::PrintRandomStatus(std::ostream& ostr) const
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Herwig7Interface::SetUserParameters()
{
  G4cout << "set user parameters of Herwig common." << G4endl
         << "nothing to be done in default."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMC::GenEvent* HepMCG4Herwig7Interface::GenerateHepMCEvent()
{

    bool isLambdaC = false;
    bool isLambda = false;
    //while(!(isLambdaC && isLambda)) {
    while(!isLambdaC) {
        m_event = m_hw->shoot();
        ThePEG::tcPVector all;
        m_event->select(std::back_inserter(all), ThePEG::SelectAll());
        std::stable_sort(all.begin(), all.end(), ThePEG::ParticleOrderNumberCmp());

        for (int i = 0; i < all.size(); ++i){
            ThePEG::tcPPtr p = all[i];
            if(p->id() == 4122) {
                isLambdaC = true;
                break;
            }
        /*    if(pyth_event[i].id() == 3122) {
                isLambda = true;
                //break;
            } */
        }

        /*if(isLambdaC && isLambda) {
            break;
        }
        else {
            isLambda = false;
            isLambdaC = false;
        }*/
    }
   HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::MEV, HepMC::Units::MM);
   ThePEG::HepMCConverter<HepMC::GenEvent>::convert(*m_event, *hepmcevt, false, ThePEG::MeV, ThePEG::millimeter);
   if(verbose>0) hepmcevt-> print();

  return hepmcevt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Herwig7Interface::Print() const
{
  G4cout << "Herwig7Interface::Print()..." << G4endl;
}

//#endif
