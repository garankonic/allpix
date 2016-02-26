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
// $Id: AllPixExtDecayerPhysics.hh 72244 2013-07-12 08:49:56Z gcosmo $
// 
/// \file eventgenerator/pythia/decayer8/include/AllPixExtDecayerPhysics.hh
/// \brief Definition of the AllPixExtDecayerPhysics class
///
/// \author I. Hrivnacova; IPN Orsay

#ifndef AllPix_EXT_DECAYER_PHYSICS_H
#define AllPix_EXT_DECAYER_PHYSICS_H

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4Decay;

/// The builder for external decayer.
///
/// The external decayer is added to all instantiated decay
/// processes
///
/// \author I. Hrivnacova; IPN Orsay

class AllPixExtDecayerPhysics: public G4VPhysicsConstructor
{
  public:
    AllPixExtDecayerPhysics(const G4String& name = "ExtDecayer");
    virtual ~AllPixExtDecayerPhysics();
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  protected:


  private:
    /// Not implemented
    AllPixExtDecayerPhysics(const AllPixExtDecayerPhysics& right);
    /// Not implemented
    AllPixExtDecayerPhysics& operator=(const AllPixExtDecayerPhysics& right);
};

#endif //AllPix_EXT_DECAYER_PHYSICS_H

