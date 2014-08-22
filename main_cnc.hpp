/** header file for main.cpp
 * **/
#ifndef __MAIN_H__
#define __MAIN_H__

#include "types.h"

#if defined (CNC) || (OMP) 
//call 2D_main from main.cpp
//template <typename T> void main_2DKriging(Input in, T context);
void main_2DKriging(Input in);
#endif//CNC||OMP

#endif //__MAIN_H__
