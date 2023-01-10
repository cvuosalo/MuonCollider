#ifndef BIBModule_h
#define BIBModule_h

/** \class BIBModule
 * 
 *  Psudo-simulation of BIB effect with given distribution of energy and position.
 * 
 *  \author H. Jia - UW-Madison
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class BIBModule: public DelphesModule
{
public:
  BIBModule();
  ~BIBModule();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItInputArray; //!
  TIterator *fItStableInputArray; //!

  const TObjArray *fInputArray; //!
  const TObjArray *fStableInputArray; //!

  TObjArray *fOutputArray; //!
  TObjArray *fStableOutputArray; //!

  String_t xHistName; //!
  String_t PositionHistName; //!
  String_t PhiHistName; //!
  String_t ThetaHistName; //!
  String_t MomentumHistName; //!
  String_t PdgEnergyHistName; //!

  String_t FileName; //!

  Int_t *fNumParticles; //!

  ClassDef(BIBModule, 1)
						  
};

#endif
