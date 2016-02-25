// @(#)geant4/pythia8:$Name$:$Id$
// Author: Andreas Morsch   27/10/2007

#ifndef PYTHIA_G4Pythia8
#define PYTHIA_G4Pythia8

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// G4Pythia8                                                                   //
//                                                                            //
// TPythia is an interface class to C++ version of Pythia 8.1                 //
// event generators, written by T.Sjostrand.                                  //
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

#include "Pythia8/Pythia.h"
#include "globals.hh"

class G4Pythia8
{
private:

protected:
   static  G4Pythia8       *fgInstance;             //! singleton instance
   Pythia8::Pythia        *fPythia;                //! The pythia8 instance
   G4int                   fNumberOfParticles;     //! Number of particles
public:
   G4Pythia8();
   G4Pythia8(const char *xmlDir);
   virtual ~G4Pythia8();
   static G4Pythia8        *Instance();
   Pythia8::Pythia        *Pythia8() {return fPythia;}

   // Interface
   virtual void            GenerateEvent();

   // Others
   void                    ReadString(const char* string) const;
   void                    ReadConfigFile(const char* string) const;
   void                    SetSeed(G4int iseed) const;
   G4bool                  Initialize(G4int idAin, G4int idBin, G4double ecms);
   G4bool                  Initialize(G4int idAin, G4int idBin, G4double eAin, G4double eBin);
   void                    ListAll() const;
   void                    ListChanged() const;
   void                    Plist(G4int id) const;
   void                    PlistAll() const;
   void                    PlistChanged() const;
   void                    PrintStatistics() const;
   void                    EventListing() const;
   G4int                   GetN() const;

};

#endif
