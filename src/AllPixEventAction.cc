/**
 *  Author John Idarraga <idarraga@cern.ch>
 */

#include "AllPixEventAction.hh"
#include "G4DigiManager.hh"
#include "G4SDManager.hh"
#include "G4PrimaryVertex.hh"
#include "AllPixMimosa26Digitizer.hh"
#include "AllPixFEI3StandardDigitizer.hh"
#include "G4RunManager.hh"
// geometry
#include "ReadGeoDescription.hh"

AllPixEventAction::AllPixEventAction(AllPixRunAction* run){

	m_run_action = run;
	m_nHC = 0;
	m_nDigitizers = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixEventAction::~AllPixEventAction()
{

}



G4String AllPixEventAction::GetNewName(G4String oldName, G4String toErase, G4String toAppend){

	string oldName_s = oldName.data();
	string toErase_s = toErase.data();
	string toAppend_s = toAppend.data();
	string fixedString;

	size_t pos = oldName_s.find_last_of(toErase) - toErase_s.size();
	if(pos == string::npos) // couldn't find it
		return oldName;


	fixedString = oldName_s.substr(0,pos);
	fixedString.append(toAppend_s);

	return G4String(fixedString);
}

void AllPixEventAction::SetDetectorDigitInput(G4double thl, G4int detId){

	m_digiPtrs[detId]->SetDetectorDigitInputs(thl);

}

void AllPixEventAction::BeginOfEventAction(const G4Event * /*evt*/)
{  

	//G4PrimaryVertex * pv = evt->GetPrimaryVertex();
    m_run_action->track_pdgid.clear();

}

void AllPixEventAction::EndOfEventAction(const G4Event * evt)
{

	G4DigiManager * fDM = G4DigiManager::GetDMpointer();
	// find digitizer module and digitize
	AllPixMimosa26Digitizer * myDM;

	G4PrimaryVertex * pv = evt->GetPrimaryVertex();

	for(G4int itr = 0 ; itr < m_nDigitizers ; itr++){
		myDM = (AllPixMimosa26Digitizer*)fDM->FindDigitizerModule( digitizerModulesNames[itr] );
		myDM->SetPrimaryVertex(pv);
		myDM->Digitize();
	}

    /*G4cout<<"!!!!!! PDG ID of each track\ntrack\tpdg\n";
    for (std::map<int,pair<int,int>>::iterator it=m_run_action->track_pdgid.begin(); it!=m_run_action->track_pdgid.end(); ++it)
    {    G4cout << it->first << "\t" << it->second.first << "\t" << it->second.second << "\n";
    }*/
    m_run_action->track_pdgid.clear();

	// digits will be retrieved at the end of the event in AllPixRun.

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

