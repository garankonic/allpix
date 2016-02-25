// @(#)root/pythia8:$Name$:$Id$
// Author: Andreas Morsch   27/10/2007

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// G4Pythia8                                                                   //
//                                                                            //
// TPythia is an interface class to C++ version of Pythia 8.1                 //
// event generators, written by T.Sjostrand.                                  //
//                                                                            //
// The user is assumed to be familiar with the Pythia package.                //
// This class includes only a basic interface to Pythia8. Because Pythia8 is  //
// also written in C++, its functions/classes can be called directly from a   //
// compiled C++ script.                                                       //
// To call Pythia functions not available in this interface a dictionary must //
// be generated.                                                              //
// see $ROOTSYS/tutorials/pythia/pythia8.C for an example of use from CINT.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
/*
*------------------------------------------------------------------------------------*
 |                                                                                    |
 |  *------------------------------------------------------------------------------*  |
 |  |                                                                              |  |
 |  |                                                                              |  |
 |  |   PPP   Y   Y  TTTTT  H   H  III    A      Welcome to the Lund Monte Carlo!  |  |
 |  |   P  P   Y Y     T    H   H   I    A A     This is PYTHIA version 8.100      |  |
 |  |   PPP     Y      T    HHHHH   I   AAAAA    Last date of change: 20 Oct 2007  |  |
 |  |   P       Y      T    H   H   I   A   A                                      |  |
 |  |   P       Y      T    H   H  III  A   A    Now is 27 Oct 2007 at 18:26:53    |  |
 |  |                                                                              |  |
 |  |   Main author: Torbjorn Sjostrand; CERN/PH, CH-1211 Geneva, Switzerland,     |  |
 |  |     and Department of Theoretical Physics, Lund University, Lund, Sweden;    |  |
 |  |     phone: + 41 - 22 - 767 82 27; e-mail: torbjorn@thep.lu.se                |  |
 |  |   Author: Stephen Mrenna; Computing Division, Simulations Group,             |  |
 |  |     Fermi National Accelerator Laboratory, MS 234, Batavia, IL 60510, USA;   |  |
 |  |     phone: + 1 - 630 - 840 - 2556; e-mail: mrenna@fnal.gov                   |  |
 |  |   Author: Peter Skands; CERN/PH, CH-1211 Geneva, Switzerland,                |  |
 |  |     and Theoretical Physics Department,                                      |  |
 |  |     Fermi National Accelerator Laboratory, MS 106, Batavia, IL 60510, USA;   |  |
 |  |     phone: + 41 - 22 - 767 24 59; e-mail: skands@fnal.gov                    |  |
 |  |                                                                              |  |
 |  |   The main program reference is the 'Brief Introduction to PYTHIA 8.1',      |  |
 |  |   T. Sjostrand, S. Mrenna and P. Skands, arXiv:0710.3820                     |  |
 |  |                                                                              |  |
 |  |   The main physics reference is the 'PYTHIA 6.4 Physics and Manual',         |  |
 |  |   T. Sjostrand, S. Mrenna and P. Skands, JHEP05 (2006) 026 [hep-ph/0603175]. |  |
 |  |                                                                              |  |
 |  |   An archive of program versions and documentation is found on the web:      |  |
 |  |   http://www.thep.lu.se/~torbjorn/Pythia.html                                |  |
 |  |                                                                              |  |
 |  |   This program is released under the GNU General Public Licence version 2.   |  |
 |  |   Please respect the MCnet Guidelines for Event Generator Authors and Users. |  |
 |  |                                                                              |  |
 |  |   Disclaimer: this program comes without any guarantees.                     |  |
 |  |   Beware of errors and use common sense when interpreting results.           |  |
 |  |                                                                              |  |
 |  |   Copyright (C) 2007 Torbjorn Sjostrand                                      |  |
 |  |                                                                              |  |
 |  |                                                                              |  |
 |  *------------------------------------------------------------------------------*  |
 |                                                                                    |
 *------------------------------------------------------------------------------------*
*/

#include "G4Pythia8.hh"

G4Pythia8*  G4Pythia8::fgInstance = 0;

////////////////////////////////////////////////////////////////////////////////
/// Constructor

G4Pythia8::G4Pythia8():
    fPythia(0),
    fNumberOfParticles(0)
{
   if (fgInstance)
      std::cout << "[WARNING] There's already an instance of G4Pythia8\n";

   fPythia    = new Pythia8::Pythia();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with an xmlDir (eg "../xmldoc"

G4Pythia8::G4Pythia8(const char *xmlDir):
    fPythia(0),
    fNumberOfParticles(0)
{
   if (fgInstance)
      std::cout << "[WARNING] There's already an instance of G4Pythia8\n";

   fPythia    = new Pythia8::Pythia(xmlDir);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

G4Pythia8::~G4Pythia8()
{
   delete fPythia;
}

////////////////////////////////////////////////////////////////////////////////
/// Return an instance of G4Pythia8

G4Pythia8* G4Pythia8::Instance()
{
   return fgInstance ? fgInstance : (fgInstance = new G4Pythia8()) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization

G4bool G4Pythia8::Initialize(G4int idAin, G4int idBin, G4double ecms)
{

   // Set arguments in Settings database.
   fPythia->settings.mode("Beams:idA",  idAin);
   fPythia->settings.mode("Beams:idB",  idBin);
   fPythia->settings.mode("Beams:frameType",  1);
   fPythia->settings.parm("Beams:eCM", ecms);

   return fPythia->init();

   //return fPythia->init(idAin, idBin, ecms);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization

G4bool G4Pythia8::Initialize(G4int idAin, G4int idBin, G4double eAin, G4double eBin)
{

   // Set arguments in Settings database.
   fPythia->settings.mode("Beams:idA",  idAin);
   fPythia->settings.mode("Beams:idB",  idBin);
   fPythia->settings.mode("Beams:frameType",  2);
   fPythia->settings.parm("Beams:eA",      eAin);
   fPythia->settings.parm("Beams:eB",      eBin);

   // Send on to common initialization.
   return fPythia->init();

   //return fPythia->init(idAin, idBin, eAin, eBin);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate the next event

void G4Pythia8::GenerateEvent()
{
   fPythia->next();
   fNumberOfParticles  = fPythia->event.size() - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization

G4int G4Pythia8::GetN() const
{
   return (fPythia->event.size() - 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Configuration

void G4Pythia8::ReadString(const char* string) const
{
   fPythia->readString(string);
}

////////////////////////////////////////////////////////////////////////////////
/// Configuration

void  G4Pythia8::ReadConfigFile(const char* string) const
{
  fPythia->readFile(string);
}
////////////////////////////////////////////////////////////////////////////////
/// Set Seed

void G4Pythia8::SetSeed(G4int iseed) const
{
    fPythia->readString("Random:setSeed = on");
    ostringstream Seed;
    Seed<<"Random:seed = "<<iseed;
    fPythia->readString(Seed.str());
}

////////////////////////////////////////////////////////////////////////////////
/// Event listing

void G4Pythia8::ListAll() const
{
   fPythia->settings.listAll();
}

////////////////////////////////////////////////////////////////////////////////
/// Event listing

void G4Pythia8::ListChanged() const
{
   fPythia->settings.listChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Event listing

void G4Pythia8::Plist(G4int id) const
{
   fPythia->particleData.list(id);
}

////////////////////////////////////////////////////////////////////////////////
/// Event listing

void G4Pythia8::PlistAll() const
{
   fPythia->particleData.listAll();
}

////////////////////////////////////////////////////////////////////////////////
/// Event listing

void G4Pythia8::PlistChanged() const
{
   fPythia->particleData.listChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Print end of run statistics

void G4Pythia8::PrintStatistics() const
{
   fPythia->stat();
}

////////////////////////////////////////////////////////////////////////////////
/// Event listing

void G4Pythia8::EventListing() const
{
   fPythia->event.list();
}

