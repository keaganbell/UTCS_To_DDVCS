#ifndef SimpleField_H
#define SimpleField_H 1

/**
  * @file
  * @brief Defines class SimpleField
*/

#include "globals.hh"
#include "G4MagneticField.hh"
 
#include <fstream>
#include <vector>

using namespace std;

class SimpleField 

#ifndef STANDALONE
: public G4MagneticField
#endif

{
public:
  SimpleField();
  ~SimpleField();
  void GetFieldValue( const  double Point[4],
		      double *Bfield ) const;
  void TestFieldMap(double xlo, double xhi, double zlo, double zhi,
		    double xstp, double zstp);
  void TestFieldMap(double xlo, double xhi,
		    double ylo, double yhi,
		    double zlo, double zhi,
		    double xstp, double ystp, double zstp);

private:
  G4double fStrength;
  G4double fXCenter, fZCenter;
  G4double fSinRotA, fCosRotA;
  bool fDebugFlag;
  /*  
  vector< vector< vector< double > > > fXField;
  vector< vector< vector< double > > > fYField;
  vector< vector< vector< double > > > fZField;
  G4int fNx,fNy,fNz; 
  G4double fMinix, fMaxix, fMiniy, fMaxiy, fMiniz, fMaxiz;
  G4double fDx, fDy, fDz;
  */
  vector< vector< double > > fRField;
  vector< vector< double > > fZField;
  G4int fNr,fNz; 
  G4double fMinir, fMaxir, fMiniz, fMaxiz;
  G4double fDr, fDz;
  
};

#endif
