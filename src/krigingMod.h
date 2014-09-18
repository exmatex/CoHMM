/** header file for krigingMod.ci
 * **/
#ifndef __KRMOD_H__
#define __KRMOD_H__
#include "krigingMod.decl.h"

class krigingChare : public CBase_krigingChare {

 public:

  krigingChare_SDAG_CODE

  /// Constructors ///
  krigingChare() {}
  krigingChare(CkMigrateMessage* msg) { }

  /// Entry Methods ///
};

class fluxOutMsg : public CMessage_fluxOutMsg {

  public:
    fluxOutput* fluxOut;

};

#endif //__KRMOD_H__
