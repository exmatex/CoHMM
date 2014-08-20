/** header file for main.C
 * **/
#ifndef __MAIN_H__
#define __MAIN_H__
#include "krigingMod.decl.h"
#include "main.decl.h"

class Main : public CBase_Main {
  CProxy_krigingChare krigingChareProxy;

 public:

  /// Constructors ///
  Main(CkArgMsg* msg);
  Main(CkMigrateMessage* msg);

  /// Entry Methods ///
  void go(Input in);

};
#endif //__MAIN_H__
