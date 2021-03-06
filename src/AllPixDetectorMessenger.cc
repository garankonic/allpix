/**
 *  Author John Idarraga <idarraga@cern.ch>
 */
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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "AllPixDetectorMessenger.hh"

#include "AllPixDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixDetectorMessenger::AllPixDetectorMessenger(
		AllPixDetectorConstruction* AllPixDet)
: m_AllPixDetector(AllPixDet)
{ 
	m_allpixDir = new G4UIdirectory("/allpix/");
	m_allpixDir->SetGuidance("UI commands MedipixG4");

	m_detDir = new G4UIdirectory("/allpix/det/");
	m_detDir->SetGuidance("detector control");

	m_ROOTDir = new G4UIdirectory("/allpix/WriteROOTFiles/");
	m_ROOTDir->SetGuidance("ROOT File output control");

	m_configDir = new G4UIdirectory("/allpix/config/");
	m_configDir->SetGuidance("allpix configuration");

	m_extrasDir = new G4UIdirectory("/allpix/extras/");
	m_extrasDir->SetGuidance("extras");

	m_beamDir = new G4UIdirectory("/allpix/beam/");
	m_beamDir->SetGuidance("beam");

#ifdef _EUTELESCOPE
	m_detDir = new G4UIdirectory("/allpix/eudet/");
	m_detDir->SetGuidance("EUDET");
#endif

	m_detIdCmd = new G4UIcmdWithAnInteger("/allpix/det/setId", this);
	m_detIdCmd->SetGuidance("Detector ID");
	m_detIdCmd->SetParameterName("ID", true);
	m_detIdCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_detPosCmd = new G4UIcmdWith3VectorAndUnit("/allpix/det/setPosition", this);
	m_detPosCmd->SetGuidance("Set position of the center of the Si wafer,");
	m_detPosCmd->SetGuidance("If you call this function N times there will be N detectors.");
	m_detPosCmd->SetParameterName("posx", "posy", "posz", false, false); // non omittable, no default
	m_detPosCmd->SetUnitCategory("Length");
	m_detPosCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_detRotCmd = new G4UIcmdWith3VectorAndUnit("/allpix/det/setRotation", this);
	m_detRotCmd->SetGuidance("Set rotation of a medipix.  If you don't call this command the medipix");
	m_detRotCmd->SetGuidance("will be placed without rotation.");
	m_detRotCmd->SetParameterName("rotx", "roty", "rotz", false, false); // non omittable, no default
	m_detRotCmd->SetUnitCategory("Angle");
	m_detRotCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	///////////
	m_UpdateCmd = new G4UIcmdWithoutParameter("/allpix/det/update",this);
	m_UpdateCmd->SetGuidance("Update geometry.");
	m_UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
	m_UpdateCmd->SetGuidance("if you changed geometrical value(s).");
	m_UpdateCmd->AvailableForStates(G4State_Idle);
	///////////

	m_HighTHLCmd = new G4UIcmdWithADoubleAndUnit("/allpix/det/setHighTHL",this);
	m_HighTHLCmd->SetGuidance("High Threshold.");
	m_HighTHLCmd->SetParameterName("highTHL", false, false);
	m_HighTHLCmd->SetUnitCategory("Energy");
	m_HighTHLCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_LowTHLCmd = new G4UIcmdWithADoubleAndUnit("/allpix/det/setLowTHL",this);
	m_LowTHLCmd->SetGuidance("Low Threshold.");
	m_LowTHLCmd->SetParameterName("lowTHL", false, false);
	m_LowTHLCmd->SetUnitCategory("Energy");
	m_LowTHLCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_AcqTimeCmd = new G4UIcmdWithADoubleAndUnit("/allpix/det/setAcqTime",this);
	m_AcqTimeCmd->SetGuidance("Aquisition time.");
	m_AcqTimeCmd->SetParameterName("acqTime", false, false);
	m_AcqTimeCmd->SetUnitCategory("Time");
	m_AcqTimeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_HVCmd = new G4UIcmdWithADoubleAndUnit("/allpix/det/setHV",this);
	m_HVCmd->SetGuidance("Bias voltage.");
	m_HVCmd->SetParameterName("HV", false, false);
	//m_HVCmd->SetUnitCategory("");
	m_HVCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_ClockCmd = new G4UIcmdWithADoubleAndUnit("/allpix/det/setClock",this);
	m_ClockCmd->SetGuidance("The clock.");
	m_ClockCmd->SetParameterName("Clock", false, false);
	m_ClockCmd->SetUnitCategory("Frequency");
	m_ClockCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_StepLengthSensor = new G4UIcmdWithADoubleAndUnit("/allpix/det/setMaxStepLengthSensor",this);
	m_StepLengthSensor->SetGuidance("User Limit. Max step length. Applies only to the sensor.");
	m_StepLengthSensor->SetParameterName("MaxStepLengthSensor", false, false);
	m_StepLengthSensor->SetUnitCategory("Distance");
	m_StepLengthSensor->AvailableForStates(G4State_PreInit, G4State_Idle);

	//////////////////////////
	// Config

	m_outputPrefix = new G4UIcmdWithAString("/allpix/config/setOutputPrefixWithPath", this);
	m_outputPrefix->SetGuidance("Set output file prefix (path can be included).  If no path is specified, the file will be written to ./");
	m_outputPrefix->SetParameterName("OutputPrefix", false);
	m_outputPrefix->SetDefaultValue("allpixoutput");
	m_outputPrefix->AvailableForStates(G4State_PreInit, G4State_Idle);

	//////////////////////////
	// extras

	// Appliance structure placed with respect to the wrapper (medipix)
	m_detAppliancePosCmd = new G4UIcmdWith3VectorAndUnit("/allpix/extras/setAppliancePosition", this);
	m_detAppliancePosCmd->SetGuidance("Set position of the detector appliance.  If you don't call this command the");
	m_detAppliancePosCmd->SetGuidance("structure will not be built.");
	m_detAppliancePosCmd->SetParameterName("posx", "posy", "posz", false, false); // non omittable, no default
	m_detAppliancePosCmd->SetUnitCategory("Length");
	m_detAppliancePosCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_wrapperEnhancementCmd = new G4UIcmdWith3VectorAndUnit("/allpix/extras/setWrapperEnhancement", this);
	m_wrapperEnhancementCmd->SetGuidance("Enhace the wrapper volume to include appliances, if needed");
	m_wrapperEnhancementCmd->SetParameterName("offsetx", "offsety", "offsetz", false, false); // non omittable, no default
	m_wrapperEnhancementCmd->SetUnitCategory("Length");
	m_wrapperEnhancementCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_AppliancesType = new G4UIcmdWithAnInteger("/allpix/extras/setApplianceType", this);
	m_AppliancesType->SetParameterName("Type", true);
	//m_ApplianceType->SetDefaultValue(0);
	m_AppliancesType->AvailableForStates(G4State_PreInit);

	// Test structure placed with respect to the world
	m_testStructPosCmd = new G4UIcmdWith3VectorAndUnit("/allpix/extras/setTestStructurePosition", this);
	m_testStructPosCmd->SetGuidance("Set position of the test structure.  If you don't call this command the");
	m_testStructPosCmd->SetGuidance("test structure will not be built.");
	m_testStructPosCmd->SetParameterName("posx", "posy", "posz", false, false); // non omittable, no default
	m_testStructPosCmd->SetUnitCategory("Length");
	m_testStructPosCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_testStructRotCmd = new G4UIcmdWith3VectorAndUnit("/allpix/extras/setTestStructureRotation", this);
	m_testStructRotCmd->SetGuidance("Set rotation of the test structure.");
	m_testStructRotCmd->SetParameterName("rotx", "roty", "rotz", false, false); // non omittable, no default
	m_testStructRotCmd->SetUnitCategory("Angle");
	m_testStructRotCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_testStructDetLinkCmd = new G4UIcmdWithAnInteger("/allpix/extras/setTestStructureDetLink", this);
	m_testStructDetLinkCmd->SetGuidance("Detector ID for related sensor");
	m_testStructDetLinkCmd->SetParameterName("ID", true);
	m_testStructDetLinkCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_testStructType = new G4UIcmdWithAnInteger("/allpix/extras/setTestStructureType", this);
	m_testStructType->SetParameterName("Type", true);
	m_testStructType->SetDefaultValue(9999);
	m_testStructType->AvailableForStates(G4State_PreInit);

	// Set material for world volume
	m_worldMaterial = new G4UIcmdWithAString("/allpix/extras/setWorldMaterial", this);
	m_worldMaterial->SetGuidance("Set material for world volume. (string: \"Air\", \"Vacuum\")");
	m_worldMaterial->SetParameterName("Material", false); // non omittable, no default
	m_worldMaterial->SetDefaultValue("Air");
	m_worldMaterial->SetCandidates("Air Vacuum");
	m_worldMaterial->AvailableForStates(G4State_PreInit, G4State_Idle);

	m_magFieldCmd = new G4UIcmdWithADoubleAndUnit("/allpix/extras/setPeakField",this);
	m_magFieldCmd->SetGuidance("Define magnetic field peak value.");
	m_magFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
	m_magFieldCmd->SetParameterName("Bz",false);
	m_magFieldCmd->SetUnitCategory("Magnetic flux density");
	m_magFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fMyXtalDirectory = new G4UIdirectory("/xtal/");
    fMyXtalDirectory->SetGuidance("Xtal setup control commands.");

    fXtalMaterialCmd = new G4UIcmdWithAString("/xtal/setMaterial",
                                              this);
    fXtalMaterialCmd->SetGuidance("Set Xtal material.");
    fXtalMaterialCmd->SetParameterName("xMat",true);
    fXtalMaterialCmd->SetDefaultValue("G4_Si");

    fXtalCurvatureRadiusCmd =
        new G4UIcmdWith3VectorAndUnit("/xtal/setCurvRadius",this);
    fXtalCurvatureRadiusCmd->SetGuidance("Set Xtal curvature radius.");
    fXtalCurvatureRadiusCmd->SetParameterName("xtalcrx",
                                              "xtalcry",
                                              "xtalcrz",
                                              true);
    fXtalCurvatureRadiusCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    fXtalCurvatureRadiusCmd->SetDefaultUnit("m");

    fXtalSizeCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setSize",this);
    fXtalSizeCmd->SetGuidance("Set Xtal size.");
    fXtalSizeCmd->SetParameterName("xtalSizeX",
                                   "xtalSizeY",
                                   "xtalSizeZ",
                                   true);
    fXtalSizeCmd->SetDefaultValue(G4ThreeVector(6.,2.,6.));
    fXtalSizeCmd->SetDefaultUnit("mm");

    fXtalAngleCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setAngle",this);
    fXtalAngleCmd->SetGuidance("Set Xtal angles in the world.");
    fXtalAngleCmd->SetParameterName("xtalAngleX",
                                    "xtalAngleY",
                                    "xtalAngleZ",
                                    true);
    fXtalAngleCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    fXtalAngleCmd->SetDefaultUnit("rad");

    fXtalCellSizeCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setCellSize",this);
    fXtalCellSizeCmd->SetGuidance("Set Xtal cell size.");
    fXtalCellSizeCmd->SetParameterName("xtalCellSizeX",
                                       "xtalCellSizeY",
                                       "xtalCellSizeZ",
                                       true);
    fXtalCellSizeCmd->SetDefaultValue(G4ThreeVector(5.431,5.431,5.431));
    fXtalCellSizeCmd->SetDefaultUnit("angstrom");

    fXtalCellAngleCmd =
        new G4UIcmdWith3VectorAndUnit("/xtal/setCellAngle",this);
    fXtalCellAngleCmd->SetGuidance("Set Xtal cell angles in the world.");
    fXtalCellAngleCmd->SetParameterName("xtalCellAngleX",
                                        "xtalCellAngleY",
                                        "xtalCellAngleZ",
                                        true);
    fXtalCellAngleCmd->SetDefaultValue(G4ThreeVector(90.,90.,90.));
    fXtalCellAngleCmd->SetDefaultUnit("deg");

    fXtalCellThermalVibration =
        new G4UIcmdWithADoubleAndUnit("/xtal/setThermVibr",this);
    fXtalCellThermalVibration->SetGuidance("Thermal vibration amplitude");
    fXtalCellThermalVibration->SetParameterName("xtalThermVibr",true);
    fXtalCellThermalVibration->SetDefaultValue(0.075);
    fXtalCellThermalVibration->SetDefaultUnit("angstrom");

    fXtalMillerCmd = new G4UIcmdWith3Vector("/xtal/setMiller",this);
    fXtalMillerCmd->SetGuidance("Set Miller Indexes");
    fXtalMillerCmd->SetParameterName("xtalMiller1",
                                     "xtalMiller2",
                                     "xtalMiller3",
                                     true);
    fXtalMillerCmd->SetDefaultValue(G4ThreeVector(2,2,0));

    fAddXtal = new G4UIcmdWithAString("/addxtal",this);
    fAddXtal->SetGuidance("Add Xtal.");
    fAddXtal->SetParameterName("addxtal",true);



#ifdef _EUTELESCOPE
	// Specific EUTelescope
	m_scint1PosCmd = new G4UIcmdWith3VectorAndUnit("/allpix/eudet/scint1Pos", this);
	m_scint2PosCmd = new G4UIcmdWith3VectorAndUnit("/allpix/eudet/scint2Pos", this);
	m_scint3PosCmd = new G4UIcmdWith3VectorAndUnit("/allpix/eudet/scint3Pos", this);
	m_scint4PosCmd = new G4UIcmdWith3VectorAndUnit("/allpix/eudet/scint4Pos", this);
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixDetectorMessenger::~AllPixDetectorMessenger()
{

	delete m_detPosCmd;
	delete m_detRotCmd;
	delete m_testStructPosCmd;
	delete m_testStructRotCmd;
	delete m_detAppliancePosCmd;
	delete m_UpdateCmd;
	delete m_worldMaterial;

	delete m_outputPrefix;

	delete m_detDir;
	delete m_allpixDir;
	delete m_ROOTDir;

    delete fXtalMaterialCmd;
    delete fXtalMillerCmd;
    delete fXtalCurvatureRadiusCmd;
    delete fXtalSizeCmd;
    delete fXtalAngleCmd;
    delete fXtalCellSizeCmd;
    delete fXtalCellAngleCmd;
    delete fAddXtal;


#ifdef _EUTELESCOPE
	delete m_scint1PosCmd;
	delete m_scint2PosCmd;
	delete m_scint3PosCmd;
	delete m_scint4PosCmd;
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 

	//G4cout << "Setting up : " << newValue << G4endl;

	if( command == m_detIdCmd )
	{
		m_AllPixDetector->SetDetectorID(
				m_detIdCmd->GetNewIntValue(newValue)
		);

	}

	if( command == m_detPosCmd )
	{
		m_AllPixDetector->SetDetectorPosition(
				m_detPosCmd->GetNew3VectorValue(newValue)
		);
	}
	if( command == m_detRotCmd )
	{
		m_AllPixDetector->SetDetectorRotation(
				m_detPosCmd->GetNew3VectorValue(newValue)
		);
	}
	if( command == m_UpdateCmd )
	{ m_AllPixDetector->UpdateGeometry(); }

	if( command == m_LowTHLCmd )
	{
		m_AllPixDetector->SetLowTHL(
				m_LowTHLCmd->GetNewDoubleValue(newValue)
		);
	}

	if( command == m_testStructPosCmd )
	{
		m_AllPixDetector->SetBuildTestStructure(true);
		m_AllPixDetector->SetTestStructurePosition(
				m_testStructPosCmd->GetNew3VectorValue(newValue)
		);
	}

	if (command == m_testStructType)
	{
		m_AllPixDetector->SetTestStructureType(
				m_testStructType->GetNewIntValue(newValue)
		);
	}

	if( command == m_testStructRotCmd )
	{
		m_AllPixDetector->SetTestStructureRotation(
				m_testStructRotCmd->GetNew3VectorValue(newValue)
		);
	}

	if( command == m_testStructDetLinkCmd )
	{
		m_AllPixDetector->SetTestStructureDetectorLink(
				m_testStructDetLinkCmd->GetNewIntValue(newValue)
		);

	}

	if( command == m_AppliancesType)
	{
		m_AllPixDetector->SetAppliancesType(
				m_AppliancesType->GetNewIntValue(newValue)
				);
	}

	if( command == m_detAppliancePosCmd )
	{
		m_AllPixDetector->SetBuildAppliances(true);
		m_AllPixDetector->SetAppliancePosition(
				m_detAppliancePosCmd->GetNew3VectorValue(newValue)
		);
	}

	if( command == m_wrapperEnhancementCmd )
	{
		m_AllPixDetector->SetWrapperEnhancement(
				m_wrapperEnhancementCmd->GetNew3VectorValue(newValue)
		);
	}

	if( command == m_worldMaterial )
	{
		m_AllPixDetector->SetWorldMaterial(
				newValue
		);
	}

	if( command == m_magFieldCmd )
	{
		G4cout << "Setting up magnetic field " << m_magFieldCmd->GetNewDoubleValue(newValue) << G4endl;
		m_AllPixDetector->SetPeakMagField(
				m_magFieldCmd->GetNewDoubleValue(newValue)
				);
	}

	if( command == m_outputPrefix ) 
	  {
	    G4cout << "Setting up output file prefix " << newValue << G4endl;
	    m_AllPixDetector->SetOutputFilePrefix( newValue );
	  }

	if( command == m_StepLengthSensor) {
		G4cout << "Setting up Max Step Length" << newValue << G4endl;
		m_AllPixDetector->SetMaxStepLengthSensor(
				m_StepLengthSensor->GetNewDoubleValue(newValue)
		);
	}
    if(command==fXtalMaterialCmd ){
        m_AllPixDetector->SetXtalMaterial(newValue);
    }
    if(command==fXtalCurvatureRadiusCmd ){
        m_AllPixDetector->SetXtalCurvatureRadius(
                        fXtalCurvatureRadiusCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalSizeCmd ){
        m_AllPixDetector->SetXtalSize(fXtalSizeCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalAngleCmd ){
        m_AllPixDetector->SetXtalAngle(fXtalAngleCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalCellSizeCmd ){
        m_AllPixDetector->SetXtalCellSize(
                        fXtalCellSizeCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalCellAngleCmd ){
        m_AllPixDetector->SetXtalCellAngle(
                        fXtalCellAngleCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalCellThermalVibration ){
        m_AllPixDetector->SetXtalThermalVibrationAmplitude(
                        fXtalCellThermalVibration->GetNewDoubleValue(newValue));
    }
    if(command==fXtalMillerCmd ){
        m_AllPixDetector->SetXtalMiller(fXtalMillerCmd->GetNew3VectorValue(newValue));
    }
    if(command==fAddXtal ){
        m_AllPixDetector->AddXtalTarget();
    }

#ifdef _EUTELESCOPE
	if( command == m_scint1PosCmd )
	{
		m_AllPixDetector->SetScintPos(
				m_scint1PosCmd->GetNew3VectorValue(newValue)
		);
	}
	if( command == m_scint2PosCmd )
	{
		m_AllPixDetector->SetScintPos(
				m_scint2PosCmd->GetNew3VectorValue(newValue)
		);
	}
	if( command == m_scint3PosCmd )
	{
		m_AllPixDetector->SetScintPos(
				m_scint3PosCmd->GetNew3VectorValue(newValue)
		);
	}
	if( command == m_scint4PosCmd )
	{
		m_AllPixDetector->SetScintPos(
				m_scint4PosCmd->GetNew3VectorValue(newValue)
		);
	}
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
