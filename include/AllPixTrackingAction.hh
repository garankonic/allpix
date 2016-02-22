#ifndef AllPixTrackingAction_h
#define AllPixTrackingAction_h 1

class G4Track;
#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include "AllPixRunAction.hh"

class AllPixTrackingAction : public G4UserTrackingAction
{
  public:
    AllPixTrackingAction(AllPixRunAction *aRun);
    virtual ~AllPixTrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

  private:
    AllPixRunAction * m_run_action;
};

#endif // ALLPIXTRACKINGACTION_HH
