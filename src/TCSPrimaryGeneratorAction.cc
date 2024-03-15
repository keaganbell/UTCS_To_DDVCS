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

#include "TCSPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
//#include <TMath.h>
#include <cmath>

//#include "TCSGen.hh"
#include "G4HEPEvtInterface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSPrimaryGeneratorAction::TCSPrimaryGeneratorAction() :
  G4VUserPrimaryGeneratorAction()
{
 // Read in default beam parameters

  ifstream file("beam_definition.txt"); // Open the file for reading.

  string line;
  istringstream iss;

  getline(file, line);  iss.str(line);
  iss >> fParticleName;
  getline(file, line);  iss.str(line);
  iss >> fPmin >> fPmax;
  getline(file, line);  iss.str(line);
  iss >> fX0 >> fY0 >> fZ0;
  getline(file, line);  iss.str(line);
  iss >> fDX >> fDY >> fDZ;
  getline(file, line);  iss.str(line);
  iss >> fThetaX >> fThetaY;
  getline(file, line);  iss.str(line);
  iss >> fdThetaX >> fdThetaY;
  getline(file, line);  iss.str(line);
  string mode_flag;
  iss >> mode_flag;

  file.close();

  fPmin *= GeV;
  fPmax *= GeV;
  fX0 *= cm;
  fY0 *= cm;
  fZ0 *= cm;
  fDX *= mm;
  fDY *= mm;
  fDZ *= mm;
  fThetaX *= degree;
  fThetaY *= degree;
  fdThetaX *= degree;
  fdThetaY *= degree;

  if (mode_flag == "tcs")
    fMode = tcs;
  else if (mode_flag == "beam")
    fMode = beam;
  else if (mode_flag == "brem")
    fMode = brem;
  else if (mode_flag == "hep_data")
    fMode = hep_data;
  else
    fMode = beam;

  G4cout << "TCSPrimaryGeneratorAction: Initial beam definition:" << G4endl;
  G4cout << "  Requested mode: " << mode_flag << endl;
  switch (fMode) {
  case beam :
    G4cout << "  Beam mode." << G4endl;
    G4cout << "  Particle " << fParticleName << G4endl;
    G4cout << "  Momentum range: " << fPmin/GeV << " -- " << fPmax/GeV << " GeV"
	   << G4endl;
    G4cout << "  Position: (" << fX0/cm << ", " << fY0/cm << ", " << fZ0/cm
	   << ") cm" << G4endl;
    G4cout << "  Beam sizes: " << fDX/mm << " x " << fDY/mm << " x " << fDZ
	   << " mm^3" << G4endl;
    G4cout << "  Beam deflection: ThetaX = " << fThetaX/degree
	   << " deg, ThetaY = " << fThetaY/degree << " deg" << G4endl;
    if (fdThetaY<0)
      G4cout << "  Beam divergence: " << fdThetaX/degree << " deg" << G4endl;
    else
      G4cout << "  Beam sampling in solid angle with dThetaX = "
	     << fdThetaX/degree << " deg,  dThetaY = "
	     << fdThetaY/degree << " deg" << G4endl;
    break;
  case brem :
    G4cout << "  Bremsstrahlung mode." << G4endl;
    G4cout << "  Bremsstrahlung photons from electron beam." << G4endl;
    G4cout << "  Bremsstrahlung photon momentum range: " << fPmin/GeV << " -- "
	   << fPmax/GeV <<" GeV" << G4endl;
    G4cout << "  Position: (" << fX0/cm << ", " << fY0/cm << ", " << fZ0/cm
	   << ") cm" << G4endl;
    G4cout << "  Beam sizes: " << fDX/mm << " x " << fDY/mm << " x " << fDZ
	   << " mm^3" << G4endl;
    G4cout << "  Beam deflection: ThetaX = " << fThetaX/degree
	   << " deg, ThetaY = " << fThetaY/degree << " deg" << G4endl;
    //   G4cout << "  Beam divergence: " << fdThetaX/degree << " deg" << G4endl;
    break;
  case tcs :
    G4cout << "  *** TCS mode: will read TCS events from input file! ***"
	   << G4endl;
    break;
  case hep_data :
    G4cout << "  *** hep_data mode: will read events from input file! ***"
	   << G4endl;
    break;
  default:
    G4cout << "  Mode not defined, assume beam mode." << G4endl;
  }
  
  if (fMode == tcs) {
    //    fTCSEntryNum = 0;
    //    fTCSPartNum = 0;
    fHEPEvt = new G4HEPEvtInterface("tcs_gen.data");
  }
  else if (fMode == hep_data) {
    fHEPEvt = new G4HEPEvtInterface("hep.data");
  }
  else if (fMode == beam) {
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle=particleTable->FindParticle(fParticleName);

    //    particleTable->DumpTable();
    //    getchar();

    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  }
  else if (fMode == brem) {
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle=particleTable->FindParticle("gamma");

    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSPrimaryGeneratorAction::~TCSPrimaryGeneratorAction()
{
if (fMode == tcs || fMode == hep_data) {
  delete fHEPEvt;
  //  fHEPEvt = 0;
}
else
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of each event.

  //  G4cout << "TCSPrimaryGeneratorAction::GeneratePrimaries: fMode ";
  //  if (fMode == tcs)
  //    G4cout << "tcs";
  //  else
  //    G4cout << "beam";
  //  G4cout << G4endl;

  G4double x = 0.;
  G4double y = 0.;
  G4double z = 0.;

  // Circular beam cross section.
  //  do {
  //    x = G4UniformRand()-.5;
  //    y = G4UniformRand()-.5;
  //  } while (x*x+y*y > .25);

  //Rectangular beam cross section.
  x = G4UniformRand()-.5;
  y = G4UniformRand()-.5;
  
  z = G4UniformRand()-.5;

  x = x*fDX + fX0;
  y = y*fDY + fY0;
  z = z*fDZ + fZ0;

  //  G4cout << "   xyz: " << x << " " << y << " " << z << G4endl;

  const G4ThreeVector BeamDir(tan(fThetaX), tan(fThetaY), 1.);

  if (fMode == tcs || fMode == hep_data) {
    fHEPEvt->SetParticlePosition(G4ThreeVector(x,y,z));
    fHEPEvt->GeneratePrimaryVertex(anEvent);
  }
  else if (fMode == beam) {
    fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));

    ////fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fPX,fPY,fPZ));
    fParticleGun->SetParticleMomentumDirection(
        GenRandomDirection(BeamDir, fdThetaX, fdThetaY));
    fParticleGun->SetParticleMomentum(fPmin + G4UniformRand()*(fPmax-fPmin));

    //Sample momentum projection on the YZ plane.
    //G4ThreeVector N=GenRandomDirection(BeamDir,fdThetaX,fdThetaY);
    //fParticleGun->SetParticleMomentumDirection(N);
    //double PYZ = fPmin + G4UniformRand()*(fPmax-fPmin);
    //double magP = PYZ/sqrt(N.getY()*N.getY()+N.getZ()*N.getZ())*
    //sqrt(N.getX()*N.getX()+N.getY()*N.getY()+N.getZ()*N.getZ());
    //fParticleGun->SetParticleMomentum(magP);
    //cout << "TCSPrimaryGeneratorAction::GeneratePrimaries: P = " <<magP<<endl;
    //cout << " N: " << N.getX() << " " << N.getY() << " " << N.getZ() << endl;
    //cout << " PYZ = " << PYZ << endl;
    //getchar();

    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  else if (fMode == brem) {
    fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));
    fParticleGun->SetParticleMomentumDirection(BeamDir);
    ////    double Pe = fEmin + G4UniformRand()*(fEmax-fEmin);
    ////    fParticleGun->SetParticleEnergy(GetBremEnergy(Ee, 10.*MeV, Ee));
    fParticleGun->SetParticleEnergy(GetBremEnergy(fPmax, fPmin, fPmax));
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double TCSPrimaryGeneratorAction::GetBremEnergy(double Ee, double Eg_min,
						double Eg_max) {

  // Sample energy of bremsstrahlung photon from primary electron of energy Ee,
  // within the range from Eg_min -- Eg_max.

  // Acceptance-rejection (von Neumann) method.
  // Sample y from f(y)=1/y*(1-y+3/4*y^2).
  // Take h(y)=1/y, such that h(y)>=f(y);
  // a) sample y from normalized 1/y.
  // b) accept y if u <f(y)/h(y), u from uniform random distribution.
  // Return y*Ee.

  double y_min = Eg_min/Ee;
  double y_max = Eg_max/Ee;
  double y;

  do {
    double u1 = G4UniformRand();
    //    y = y_min*TMath::Power(y_max/y_min, u1);
    y = y_min*pow(y_max/y_min, u1);                      //y from 1/y distrib.
    double prob = 1.-y+3./4.*y*y;                        //f(y)/(4/3*1/y)
    if (G4UniformRand() < prob) break;
  } while (true);

  //  cout << "Brem. Energy = " << Ee*y/MeV << endl;
  //  getchar();

  return Ee*y;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector TCSPrimaryGeneratorAction::GenRandomDirection(
	    const G4ThreeVector Axis, const double Theta, const double Phi) {

  //Generate random vector around Axis. Treat Theta as maximum divergence angle
  //if Phi<0. Otherwise Theta is max. deflection angle in XZ plane, Phi is
  //max. defelection angle in YZ plane.

  //  cout << "TCSPrimaryGeneratorAction::GenRandomDirection:" << endl;
  //  cout << " Axis: " << Axis(0) << " " << Axis(1) << " " << Axis(2) << endl;
  // if (Phi>=0.)
  //  cout << " divegence theta = " << Theta/degree << " deg" << endl;
  // else
  //  cout << " sampling in thetaX = " << Theta/degree << " deg,  thetaY = "
  //  << Phi/degree << " deg" << endl;

  //Random vector around Z axis.
  G4ThreeVector vrand(0.,0.,1.);
  if (Phi < 0.) {
    double costheta = 1. - G4UniformRand()*(1.-cos(Theta));
    double phi = 2.*CLHEP::pi*G4UniformRand();
    vrand.setMag(1.);
    vrand.setTheta(acos(costheta));
    vrand.setPhi(phi);
    //  cout << " cos(theta) = " << costheta << "  theta = "
    //       << acos(costheta)/degree
    //       << " deg" << "  phi = " << phi/degree << " deg" << endl;
  }
  else {
    double thetaX = (2.*G4UniformRand()-1.)*Theta;
    double thetaY = (2.*G4UniformRand()-1.)*Phi;
    vrand.setX(tan(thetaX));
    vrand.setY(tan(thetaY));
    vrand.setZ(1.);
    //    cout << " thetaX = " << thetaX/degree << " deg,  thetaY = "
    //	 << thetaY/degree << " deg" << endl;
  }

  //  cout << " vrand at Z axis: " << vrand(0) << " " << vrand(1) << " "
  //       << vrand(2) << endl;

  //Rotate around axis perpendicular to Z axis and axis vector.
  G4ThreeVector rot_axis = CLHEP::HepZHat.cross(Axis);
  double rot_angle = CLHEP::HepZHat.angle(Axis);
  vrand.rotate(rot_angle, rot_axis);

  //  cout << " Rotation axis: " << rot_axis(0) << " " << rot_axis(1) << " "
  //       << rot_axis(2) << endl;
  //  cout << " Rotation angle = " << rot_angle/degree << " deg" << endl;
  //  cout << " vrand at Axis:   " << vrand(0) << " " << vrand(1) << " "
  //       << vrand(2) << endl;
  //  getchar();

  return vrand;
}
