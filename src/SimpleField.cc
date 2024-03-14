/**
  * @file
  * @brief implements class SimpleField
*/

// GN 2015 from above model

#include "SimpleField.hh"
#include "G4SystemOfUnits.hh"
//#include <TMath.h>
#include <cmath>

//-----------------------------------------------------------------------------

SimpleField::SimpleField()
{
  // Read field param-s

  ifstream file2( "field_param.txt" ); // Open the file for reading.

  string line;
  istringstream iss;

  string filename;
  G4double rot_angle;

  getline(file2, line);  iss.str(line);
  iss >> fStrength;
  getline(file2, line);  iss.str(line);
  iss >> fXCenter;
  getline(file2, line);  iss.str(line);
  iss >> fZCenter;
  getline(file2, line);  iss.str(line);
  iss >> rot_angle;
  getline(file2, line);  iss.str(line);
  iss >> fNr;
  getline(file2, line);  iss.str(line);
  iss >> fNz;
  getline(file2, line);  iss.str(line);
  iss >> filename;
  getline(file2, line);  iss.str(line);
  iss >> fDebugFlag;

  file2.close();

  rot_angle *= degree;
  fXCenter *= cm;
  fZCenter *= cm;

  G4cout << "Field parameters:" << G4endl;
  G4cout << "fStrength = "<< fStrength << G4endl;
  G4cout << "fXCenter     = " << fXCenter/cm << " cm" << G4endl;
  G4cout << "fZCenter  = " << fZCenter/cm << " cm" << G4endl;
  G4cout << "rot_angle = " << rot_angle/degree << " deg." << G4endl;
  G4cout << "fNr = " << fNr << G4endl;
  G4cout << "fNz = " << fNz << G4endl;
  G4cout << "filename  = " << filename << G4endl;
  G4cout << "fDebugFlag= " << fDebugFlag << G4endl;
  G4cout << G4endl;

  fSinRotA = sin(rot_angle);
  fCosRotA = cos(rot_angle);

  G4cout << "      Field center at (" << fXCenter/cm << ", " << fZCenter/cm
	 << ") cm" << endl;
  G4cout << "      sin & cos of field rotation: " << fSinRotA << " " << fCosRotA
	 << endl;

  filename = "field/" + filename;

  G4cout << "\n-----------------------------------------------------------"
	 << "\n      UVA Target Magnetic Field"
	 << "\n-----------------------------------------------------------";
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... "
	 << endl; 

  ifstream file( filename.c_str() ); // Open the file for reading.

  // Set up storage space for table
  fRField.resize( fNr );
  fZField.resize( fNr );
  int ir, iz;
  for (ir=0; ir<fNr; ir++) {
    fRField[ir].resize(fNz);
    fZField[ir].resize(fNz);
  }

  // Read in the data

  for (int i=0; i<5; i++) {
    string tmps;
    file >> tmps;
    G4cout << tmps << " ";
  }
  G4cout << G4endl;

  double zval,rval,br,bz,btot;
  for (ir=0; ir<fNr; ir++) {
    for (iz=0; iz<fNz; iz++) {

      file >> zval >> rval >> bz >> br >> btot;
      if ( (iz<5 && ir==0) || (iz>fNz-5 && ir==fNr-1) ) {
	G4cout << zval << " "<< rval << " " << bz << G4endl;
      }

      if ( ir==0 && iz==0 ) {
	fMinir = rval * cm;
	fMiniz = zval * cm;
      }

      fRField[ir][iz] = br*tesla;
      fZField[ir][iz] = bz*tesla;

    }
  }
  fMaxir = rval * cm;
  fMaxiz = zval * cm;
  G4cout << "\n ---> ... done reading " << G4endl;

  //  if (fDebugFlag) {
  //    cout << "SimpleField: Field Map:" << endl;
  //    for (ir=0; ir<fNr; ir++) {
  //      for (iz=0; iz<fNz; iz++) {
  //	cout << ir << " " << iz << ": "
  //	     << fRField[ir][iz] << " " << fZField[ir][iz] << endl;
  //      }
  //    }
  //  }

  G4cout << " ---> assumed the order:  r, z, Br, Bz "
	 << "\n ---> Min values r,z: " 
	 << fMinir/cm << " " << fMiniz/cm << " cm "
	 << "\n ---> Max values r,z: " 
	 << fMaxir/cm << " " << fMaxiz/cm << " cm " << endl;

  fDr = fMaxir - fMinir;
  fDz = fMaxiz - fMiniz;
  G4cout << "\n ---> Dif values r,z (range): " 
	 << fDr/cm << " cm in R, " << fDz/cm << " cm in Z"
	 << "\n-----------------------------------------------------------"
	 << endl;

  file.close();

  // Test field map.
  //  TestFieldMap(fXCenter-5*cm, fXCenter+5*cm, fZCenter-5*cm, fZCenter+5*cm,
  //	       5.*cm, 5.*cm);
  TestFieldMap(-5*cm, +5*cm, -5*cm, +5*cm, -5*cm, 5*cm,
      	       5.*cm, 5.*cm, 5.*cm);

  G4cout << "\n-----------------------------------------------------------"
	 << "\n  Done Initializing Simple Magnetic Field! "
	 << "\n-----------------------------------------------------------"
	 <<G4endl;;
  return;
  //  getchar();
}

//-----------------------------------------------------------------------------
 
SimpleField::~SimpleField()
{;}

//-----------------------------------------------------------------------------

void SimpleField::GetFieldValue(const double point[4],double *Bfield) const
{
  ////
  ////  G4int eID = 0;
  ////  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  ////  if(evt) eID = evt->GetEventID();
  ////cout << "SimpleField: eID = " << eID << endl;
  ////  if (eID >= 254) fDebugFlag = 1;
  ////

// Interpolate field MAP (a la TabulatedField3D!)
  Bfield[0] = 0.0;
  Bfield[1] = 0.0;
  Bfield[2] = 0.0;
  Bfield[3] = 0.0;
  Bfield[4] = 0.0;
  Bfield[5] = 0.0;
  
  double x = point[0];
  double y = point[1];
  double z = point[2]; 

  if(fDebugFlag)
  cout << "SimpleField::GetFieldValue: x,y,z [cm]: " << x/cm << " " << y/cm
       << " " << z/cm << endl;

  //From G4 to TOSCA system.

  x -= fXCenter;                   //move to field center
  z -= fZCenter;
  if(fDebugFlag)
    cout << " At center x,y,z [cm]: " << x/cm << " " << y/cm << " "
	 << z/cm << endl;

  G4double tmp;                    //transform to the rotated field system
  tmp = x*fCosRotA + z*fSinRotA;
  z   = z*fCosRotA - x*fSinRotA;
  x   = tmp;
  if(fDebugFlag)
    cout << " Rotated   x,y,z [cm]: " << x/cm << " " << y/cm << " " << z/cm
	 << endl;

  //  bool flag = abs(z)<=500 && abs(y)>=300 && abs(x)<=150;   //mm
  bool flag = 1;

  // Tabulated z-s are > 0, and Bz is simmetric wrt Z
  // while Br is asimmetric wrt Z!

  bool r_flip = z<0.;  // Flip Br if z<0

  //  z = TMath::Abs(z);
  z = abs(z);

  //  double r = TMath::Sqrt(x*x + y*y);
  double r = sqrt(x*x + y*y);

  if(fDebugFlag && flag) cout << " Radius R = " << r/cm << " cm" << endl;
  
  // Check that the point is within the defined region 

  if   (r>=fMinir && r<=fMaxir &&
	z>=fMiniz && z<=fMaxiz  ) {

    // Position of given point within region, normalized to the range [0,1]
    double rfraction = (r - fMinir) / fDr;
    double zfraction = (z - fMiniz) / fDz;
    if(fDebugFlag && flag)
      cout << " Fractions of r,z: " << rfraction << " " << zfraction << endl;

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double rdindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double rlocal = ( std::modf(rfraction*(fNr-1), &rdindex));
    double zlocal = ( std::modf(zfraction*(fNz-1), &zdindex));

    // One step back if at the end of range.
    if (rfraction == 1. && rlocal == 0.) {rdindex -= 1.; rlocal = 1.;}
    if (zfraction == 1. && zlocal == 0.) {zdindex -= 1.; zlocal = 1.;}

    if(fDebugFlag && flag) {
      cout << " Local r,z: " << rlocal << " " << zlocal << endl;
      cout << " DIndices: r,z: " << rdindex << " " << zdindex << endl;
    }

    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int rindex = static_cast<int>(rdindex);
    int zindex = static_cast<int>(zdindex);
    if(fDebugFlag && flag)
      cout << " Indices: r,z: " << rindex << " " << zindex << endl;

    if(fDebugFlag && flag) {
      cout << "fRFields [T]:" << endl;
      cout << fRField[rindex  ][zindex  ]/tesla << " "
	   << rindex << " " << zindex << endl;
      cout << fRField[rindex  ][zindex+1]/tesla << " "
	   << rindex << " " << zindex+1 << endl;
      cout << fRField[rindex+1][zindex  ]/tesla << " "
	   << rindex+1 << " " << zindex << endl;
      cout << fRField[rindex+1][zindex+1]/tesla << " "
	   << rindex+1 << " " << zindex+1 << endl;
      cout << "fZFields [T]:" << endl;
      cout << fZField[rindex  ][zindex  ]/tesla << " "
	   << rindex << " " << zindex << endl;
      cout << fZField[rindex  ][zindex+1]/tesla << " "
	   << rindex << " " << zindex+1 << endl;
      cout << fZField[rindex+1][zindex  ]/tesla << " "
	   << rindex+1 << " " << zindex << endl;
      cout << fZField[rindex+1][zindex+1]/tesla << " "
	   << rindex+1 << " " << zindex+1 << endl;
    }

    // Interpolated field

    double Br =
      fRField[rindex  ][zindex  ] * (1-rlocal) * (1-zlocal) +
      fRField[rindex  ][zindex+1] * (1-rlocal) *    zlocal  +
      fRField[rindex+1][zindex  ] *    rlocal  * (1-zlocal) +
      fRField[rindex+1][zindex+1] *    rlocal  *    zlocal  ;
    Br *= fStrength;

    if (r_flip) Br = -Br;    // flip if Z<0

    double Bz =
      fZField[rindex  ][zindex  ] * (1-rlocal) * (1-zlocal) +
      fZField[rindex  ][zindex+1] * (1-rlocal) *    zlocal  +
      fZField[rindex+1][zindex  ] *    rlocal  * (1-zlocal) +
      fZField[rindex+1][zindex+1] *    rlocal  *    zlocal  ;
    Bz *= fStrength;

    if(fDebugFlag && flag)
      cout << " Bz = " << Bz/tesla << " T,   Br = " << Br/tesla << " T" << endl;

    if (r != 0.) {
      Bfield[0] = Br * x/r;
      Bfield[1] = Br * y/r;
    }
    // Otherwise the tabulated Br==0.

    Bfield[2] = Bz;

  }

  if(fDebugFlag && flag)
    cout << " Bfields [T]: " << Bfield[0]/tesla << " " << Bfield[1]/tesla
	 << " " << Bfield[2]/tesla << endl;

  //From TOSCA to G4 system.

  tmp       = Bfield[0]*fCosRotA - Bfield[2]*fSinRotA;
  Bfield[2] = Bfield[2]*fCosRotA + Bfield[0]*fSinRotA;
  Bfield[0] = tmp;
  if(fDebugFlag && flag)
    cout << " Rotated B-s [T]: " << Bfield[0]/tesla << " " << Bfield[1]/tesla
	 << " " << Bfield[2]/tesla << endl;

  //  G4cout << "Strength = " << fStrength << " B(Gauss) = "
  //	 << Bfield[0]/gauss << " " << Bfield[1]/gauss << " " << Bfield[2]/gauss
  //	 << "  (x,z) = (" << x/cm << ", " << z/cm << ") cm" << endl;

  //  G4double B =
  ////TMath::Sqrt(Bfield[0]*Bfield[0]+Bfield[1]*Bfield[1]+Bfield[2]*Bfield[2]);
  //sqrt(Bfield[0]*Bfield[0]+Bfield[1]*Bfield[1]+Bfield[2]*Bfield[2]);
  //  if (B!=0)
  //G4cout << "Field: " << B/gauss << "  " << point[0]/cm << " " << point[1]/cm
  //	   << " " << point[2]/cm << endl;
  //G4cout << "Field: " << B/gauss << "  " << point[0]/cm << " " << point[1]/cm
  //<< " " << point[2]/cm << "  " << x/cm << " " << y/cm << " " << z/cm
  //	   << endl;
}

//-----------------------------------------------------------------------------

void SimpleField::TestFieldMap(double xlo, double xhi, double zlo, double zhi,
			       double xstp, double zstp)
{
  //  G4cout << "Field Map from SimpleField::TestFieldMap:" << G4endl;
  //  cout << "  xlo = " << xlo << " xhi = " << xhi << " zlo = " << zlo << " zhi = "
  //       << zhi << " xstp = " << xstp << " zstp = " << zstp << endl;

  for (G4double x=xlo; x<=xhi; x+=xstp)
    for (G4double z=zlo; z<=zhi; z+=zstp) {
      double p[] = {x, 0., z, 0.};
      double B[6] = {};
      GetFieldValue(p, B);
      ////      G4double modB = TMath::Sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
      //      G4double modB = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
      //      G4cout << "Field: " << modB/gauss << "  " << p[0]/cm << " "
      //      	     << p[1]/cm << " " << p[2]/cm << endl;
      G4cout << "Test Field [T]: " << B[0]/tesla << "  " << B[1]/tesla << "  "
	     << B[2]/tesla << "  " << x/cm << " " << z/cm << G4endl;
      //	     << B[2]/gauss << "  " << p[0]/cm << " " << p[1]/cm << " "
      //	     << p[2]/cm << G4endl;
    }

}

//-----------------------------------------------------------------------------

void SimpleField::TestFieldMap(double xlo, double xhi,
			       double ylo, double yhi,
			       double zlo, double zhi,
			       double xstp, double ystp, double zstp)
{
  // Test field in the magnet system.

  for (G4double x=xlo; x<=xhi; x+=xstp)
    for (G4double y=ylo; y<=yhi; y+=ystp)
      for (G4double z=zlo; z<=zhi; z+=zstp) {
	double p[] = {x, y, z, 0.};
	double B[6] = {};
	GetFieldValue(p, B);
	G4cout << "Test Field [T]: " << B[0]/tesla << "  " << B[1]/tesla << "  "
	       << B[2]/tesla << "  " << x/cm << " " << y/cm << " " << z/cm
	       << G4endl;
    }

}

//-----------------------------------------------------------------------------
/*
void SimpleField::TestFieldMap(double xlo, double xhi,
			       double ylo, double yhi,
			       double zlo, double zhi,
			       double xstp, double ystp, double zstp)
{
  // Test field in the magnet system.
  // x, y, z in the magnet system.

  for (G4double x=xlo; x<=xhi; x+=xstp)
    for (G4double y=ylo; y<=yhi; y+=ystp)
      for (G4double z=zlo; z<=zhi; z+=zstp) {
	double xlab = x*fCosRotA - z*fSinRotA + fXCenter;
	double ylab = y;
	double zlab = z*fCosRotA + x*fSinRotA + fZCenter;
	double p[] = {xlab, ylab, zlab, 0.};
	double B[6] = {};
	GetFieldValue(p, B);
	G4cout << "Test Field [T]: " << B[0]/tesla << "  " << B[1]/tesla << "  "
	       << B[2]/tesla << "  " << x/cm << " " << y/cm << " " << z/cm
	       << G4endl;
    }

}
*/
