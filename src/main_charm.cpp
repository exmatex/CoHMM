/** file containing the charm main function
 * calls 2D Kriging main
 * **/
#include "main_charm.hpp"
#include "input.hpp"
#include "2DKriging.hpp"
#include "krigingMod.decl.h"

#include <cstring>

/* readonly */ CProxy_Main mainProxy;

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg) {

// Display some info about this execution
// for the user.
  CkPrintf("**************************************************\n");
  CkPrintf("**************************************************\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**   Running \"Charm 2D Kriging %d processors    **\n",
            CkNumPes());
  CkPrintf("**                                              **\n");
  CkPrintf("**                                              **\n");
  CkPrintf("**************************************************\n");
  CkPrintf("**************************************************\n");

// Set the mainProxy readonly to point to a
// proxy for the Main chare object (this
// chare object).
  mainProxy = thisProxy;

  //get input values from json file
  Input in;
  char input_file[1024];
  strcpy(input_file, msg->argv[1]);
  parse_input((string)input_file, &in);
  
  if(msg->argc >2){ 
    char host[1024];
    strcpy(host, msg->argv[2]);
    in.head_node = string(host); 
    CkPrintf("Used cmd line redis host: %s\n", in.head_node.c_str());
  }
    
  //setup char array size
  CkArrayOptions opts(in.dim_x * in.dim_y);
  //Setup round robin map
  CProxy_RRMap rrMap = CProxy_RRMap::ckNew();
  opts.setMap(rrMap);

  //create chare array
  krigingChareProxy = CProxy_krigingChare::ckNew(opts);

  //start main program
  mainProxy.go(in);
  //print exit msg after sucessful run
  CkPrintf("Exiting charm++\n");
}

// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
Main::Main(CkMigrateMessage* msg) { }

/**Charm 2DKriging main
 *
 */
void main_2DKriging(Input in, CProxy_krigingChare krigingChareProxy);


// executes main of 2D Kriging
// and exits charm after sucessful run
void Main::go(Input in) {
    main_2DKriging(in, krigingChareProxy); 
    CkExit();
}

#include "main.def.h" 
