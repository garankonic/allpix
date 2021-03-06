/**
 *  Author John Idarraga <idarraga@cern.ch>
 */

#ifndef AllPixFEI3StandardDigitizer_h
#define AllPixFEI3StandardDigitizer_h 1

// interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixFEI3StandardDigit.hh"
#include "AllPixTrackerHit.hh"
#include "G4PrimaryVertex.hh"
#include "AllPixGeoDsc.hh"
#include "TString.h"
#include "TH2D.h"
#include <map>
#include <vector>

using namespace std;

/**
 *  A simple Digitizer implementation
 */

class AllPixFEI3StandardDigitizer : public AllPixDigitizerInterface {

public:
  AllPixFEI3StandardDigitizer(G4String, G4String, G4String);
  virtual ~AllPixFEI3StandardDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double);

private:
  digitInput m_digitIn;
  AllPixFEI3StandardDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction


  TH2D *hEx;
  TH2D *hEy;
  TH2D *hEz;



  G4double MyErf(G4double x);
  G4double ComputeDriftTimeUniformField(AllPixTrackerHit *hit);
  G4double ComputeDiffusionRMS(G4double tDrift);
  G4double IntegrateGaussian(G4double xhit,G4double yhit,G4double Sigma, G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy );
  G4double ApplyTrapping(G4double tDrift, G4double Energy);
  G4double GetElectricFieldNorm(G4double x, G4double y=0, G4double z=0);
  G4double MobilityElectron(G4double x, G4double y=0, G4double z=0);
  G4double MobilityHole(G4double x, G4double y=0, G4double z=0);

  void ComputeElectricField(G4double x, G4double y=0, G4double z=0);
  vector<G4double> ComputeDriftTimeFullField(G4double x, G4double y, G4double z);
  vector<G4double> RKF5Integration(G4double x, G4double y, G4double z,G4double dt);
  G4double ComputeSubHitContribution(G4double x, G4double y, G4double z,G4double Energy);
  G4double SetDt(G4double Dt,G4double ErreurMoy);
  G4int EnergyToTOT(G4double Energy, G4double threshold);
  G4double SlimEdgeEffect(G4int nX,G4double xpos,G4double eHit);
  G4bool isSlimEdge(G4int nX, G4int nY);
  void Efield1D(G4double z);
  void Efield2D(G4double x,G4double y,G4double z);

  G4double elec;

  ///////////////////////////////////////////////////////
  //	Some constant for the charge sharing computation
  G4double mobility ;
  G4double electricFieldX;
  G4double electricFieldY;
  G4double electricFieldZ;
  G4double depletionVoltage;
  G4double depletedDepth;
  G4double resistivity;
  G4bool   bulkType;
  G4double b;
  G4double c;
  G4double Neff;

  G4double detectorThickness;
  G4double biasVoltage;
  G4double fluence;
  G4double trappingTime;
  G4double Beta_electrons;
  G4double chipNoise;
  ///////////////////////////////////////////////////////

  //Numerical accuracy of RKF5
  G4double tlow;
  G4double tup;
  G4double dtIni;
  G4double Target;

  G4double Temperature;

  // Silicon electron and hole transport constants
  G4double Default_Electron_Mobility;// Electron mobility (cm2/Vs)
  G4double Default_Hole_Mobility;// Hole mobility (cm2/Vs

  G4double Default_Electron_D; // Electron mobility (cm2/s)
  G4double Default_Hole_D;// Hole mobility (cm2/s

  //mobility dependence on electric field
  G4double Electron_AlphaField ; //[um/ns]
  G4double Electron_ThetaField ;
  G4double Electron_TempNominal  ; // [K]

  G4double Electron_Beta ;
  G4double Electron_Saturation_Velocity;

  G4double Hole_AlphaField ; //[um/ns]
  G4double Hole_ThetaField ;
  G4double Hole_TempNominal; // [K]

  G4double Hole_Beta  ;
  G4double Hole_Saturation_Velocity;

  //Some physics process switches
  G4bool doTrapping;
  G4bool doFullField;
  G4bool doSlimEdge;

  G4double epsilon;
  G4double echarge;
  G4int precision;


  //////////////////////////////////////////////////////
  // Geometry Related constants
  G4double pitchX ;
  G4double pitchY ;
  G4int nPixX;
  G4int nPixY;
  G4double chargeSharingConstant;
  G4double GRShift;
  G4int FEIX;
  G4int Sensor;

  //Tuning of the chip and counters characteristics
  G4int MipTOT;
  G4int CounterDepth;
  G4int MipCharge;
  G4double Lv1Unit;

  // Digitizer precision related parameters
  //G4double maxIntegration = 5; //

   G4double sigma;



};

#endif
