/** parse input file for 2D macro solver, inlcudes 
 * boost json file parsing
 * **/

#include <iostream>
#ifdef CHARM
#include "main_charm.hpp"
#include "main.decl.h"
#endif
#include "types.h"

#include <string>
#include <cstring>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#ifndef CHARM
#define CkPrintf printf
#endif

void parse_input(string input_file, Input *in, App *CoMD)
{
  std::ifstream file(input_file.c_str());
  boost::property_tree::ptree pt;

  if ( file )
  {
     std::stringstream buffer;
     buffer << file.rdbuf();

     file.close();
     try
     {
         boost::property_tree::read_json(buffer, pt);
     }
     catch (std::exception const& e)
     {
        std::cerr << "Something went wrong in read_json" << std::endl;
        std::cerr << e.what() << std::endl;
	exit(1);
     }
  } else {
      std::cerr << "Could not open json file " << input_file << std::endl;
      exit(1);
  }
  try
  {
    BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, pt.get_child("parameter.macro_solver"))
    {
        if(v.second.get<std::string>("id") == "dim x" ){
            in->dim_x = v.second.get<int>("value");
            CkPrintf("set dim x:                %d\n", in->dim_x);
        }
        if(v.second.get<std::string>("id") == "dim y" ){
            in->dim_y = v.second.get<int>("value");
            CkPrintf("set dim y:                %d\n", in->dim_y);
        }
        if(v.second.get<std::string>("id") == "integration steps" ){
            in->int_steps = v.second.get<int>("value");
            CkPrintf("set int steps             %d\n", in->int_steps);
        }
        if(v.second.get<std::string>("id") == "redis database" ){
            in->redis_db = v.second.get<int>("value");
            CkPrintf("set redis database:       %i\n", in->redis_db);
        }
        if(v.second.get<std::string>("id") == "flush database" ){
            in->flush_db = v.second.get<int>("value");
            CkPrintf("set flush database:       %i\n", in->flush_db);
        }
        if(v.second.get<std::string>("id") == "redis headnode" ){
            in->head_node = v.second.get<std::string>("value");
            CkPrintf("set redis headnode:       %s\n", in->head_node.c_str());
        }
        if(v.second.get<std::string>("id") == "database threshold" ){
            in->db_threshold = v.second.get<double>("value");
            CkPrintf("set db threshold:         %lf\n", in->db_threshold);
        }
        if(v.second.get<std::string>("id") == "enable kriging" ){
            in->kriging = v.second.get<int>("value");
            CkPrintf("set enable kriging:       %i\n", in->kriging);
        }
        if(v.second.get<std::string>("id") == "kriging err. threshold" ){
            in->kr_threshold = v.second.get<double>("value");
            CkPrintf("set kr threshold:         %lf\n", in->kr_threshold);
        }
        if(v.second.get<std::string>("id") == "enable kriging database" ){
            in->kriging_db = v.second.get<int>("value");
            CkPrintf("set enable kriging db:    %i\n", in->kriging_db);
        }
        if(v.second.get<std::string>("id") == "apply Gaussian noise" ){
            in->gauss_noise = v.second.get<int>("value");
            CkPrintf("set apply G. noise:       %i\n", in->gauss_noise);
        }
        if(v.second.get<std::string>("id") == "Gaussian noise" ){
            in->noise = v.second.get<double>("value");
            CkPrintf("set Gaussian noise:       %lf\n", in->noise);
        }
        if(v.second.get<std::string>("id") == "dx" ){
            in->dx = v.second.get<double>("value");
            CkPrintf("set dx:                   %lf\n", in->dx);
        }
        if(v.second.get<std::string>("id") == "dy" ){
            in->dy = v.second.get<double>("value");
            CkPrintf("set dy:                   %lf\n", in->dy);
        }
        if(v.second.get<std::string>("id") == "dt_x" ){
            in->dt_x = v.second.get<double>("value");
            CkPrintf("set dt_x:                 %lf\n", in->dt_x);
        }
        if(v.second.get<std::string>("id") == "dt_y" ){
            in->dt_y = v.second.get<double>("value");
            CkPrintf("set dt_y:                 %lf\n", in->dt_y);
        }
        if(v.second.get<std::string>("id") == "grad threshold" ){
            in->grad_threshold = v.second.get<double>("value");
            CkPrintf("set gradient threshold:   %lf\n", in->grad_threshold);
        }
        if(v.second.get<std::string>("id") == "test problem" ){
            in->test_problem = v.second.get<int>("value");
            CkPrintf("set test problem:         %i\n", in->test_problem);
        }
        if(v.second.get<std::string>("id") == "fault tolerance" ){
            in->fault_tolerance = v.second.get<int>("value");
            CkPrintf("set fault tolerance:      %i\n", in->fault_tolerance);
        }
    }
  }
  catch (std::exception const& e)
  {
      std::cerr << "Something went wrong in parsing section 'parameter.macro_solver'" << std::endl;
      std::cerr << e.what() << std::endl;
#ifdef CHARM
      CkExit();
#else
      exit(1);
#endif
  }
  try
  {
    BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, pt.get_child("parameter.mini_app"))
    {
        if(v.second.get<std::string>("id") == "CoMD on" ){
            CoMD->comd_on = v.second.get<int>("value");
            CkPrintf("set CoMD eam on:          %i\n", CoMD->comd_on);
        }
        if(v.second.get<std::string>("id") == "pot name" ){
            CoMD->pot_name = v.second.get<std::string>("value");
            CkPrintf("set CoMD pot name:        %s\n", CoMD->pot_name.c_str());
        }
        if(v.second.get<std::string>("id") == "pot type" ){
            CoMD->pot_type = v.second.get<std::string>("value");
            CkPrintf("set CoMD pot type:        %s\n", CoMD->pot_type.c_str());
        }
        if(v.second.get<std::string>("id") == "eam on" ){
            CoMD->eam_on = v.second.get<int>("value");
            CkPrintf("set CoMD eam on:          %i\n", CoMD->eam_on);
        }
        if(v.second.get<std::string>("id") == "dim x" ){
            CoMD->dim_x = v.second.get<int>("value");
            CkPrintf("set CoMD dim x:           %i\n", CoMD->dim_x);
        }
        if(v.second.get<std::string>("id") == "dim y" ){
            CoMD->dim_y = v.second.get<int>("value");
            CkPrintf("set CoMD dim y:           %i\n", CoMD->dim_y);
        }
        if(v.second.get<std::string>("id") == "dim z" ){
            CoMD->dim_z = v.second.get<int>("value");
            CkPrintf("set CoMD dim z:           %i\n", CoMD->dim_z);
        }
        if(v.second.get<std::string>("id") == "integration steps" ){
            CoMD->integration_steps = v.second.get<int>("value");
            CkPrintf("set CoMD integr. steps:   %i\n", CoMD->integration_steps);
        }
        if(v.second.get<std::string>("id") == "lattice spacing" ){
            CoMD->lattice_spacing = v.second.get<double>("value");
            CkPrintf("set CoMD latt. spacing:   %lf\n", CoMD->lattice_spacing);
        }
        if(v.second.get<std::string>("id") == "dt" ){
            CoMD->dt = v.second.get<double>("value");
            CkPrintf("set CoMD dt:              %lf\n", CoMD->dt);
        }
        CkPrintf("Using input file values for CoMD not fully implemented yet!!!\n");
        //FIXME include all other values and hand App struct to flux fkt
    }
  }
  catch (std::exception const& e)
  {
      std::cerr << "Something went wrong in parsing section 'parameter.mini_app'" << std::endl;
      std::cerr << e.what() << std::endl;
#ifdef CHARM
      CkExit();
#else
      exit(1);
#endif
  }
}
