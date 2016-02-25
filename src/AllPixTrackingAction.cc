#include "AllPixTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4RunManager.hh"
#include "map"
AllPixTrackingAction::AllPixTrackingAction(AllPixRunAction* aRun)
{
    m_run_action = aRun;
}

AllPixTrackingAction::~AllPixTrackingAction()
{;}

void AllPixTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    if(aTrack->GetParticleDefinition()->GetPDGEncoding()!=11) {
        //G4cout<<"Start tracking of "<<aTrack->GetParticleDefinition()->GetPDGEncoding()<<G4endl;
        //G4ParticleTable::GetParticleTable()->FindParticle(aTrack->GetParticleDefinition()->GetPDGEncoding())->DumpTable();

        m_run_action->track_pdgid.insert(make_pair(aTrack->GetTrackID(),make_pair(aTrack->GetParticleDefinition()->GetPDGEncoding(),aTrack->GetParentID())));
        //m_run_action->track_pdgid.insert(std::pair<int,int>(aTrack->GetTrackID(),aTrack->GetParticleDefinition()->GetPDGEncoding()));
    }
}

void AllPixTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{;}
