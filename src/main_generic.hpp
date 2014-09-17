/** header file for main.cpp
 * **/
#ifndef __MAIN_H__
#define __MAIN_H__

#include "types.h"

#if defined (CNC) || (OMP) || (CIRCLE)
//call 2D_main from main.cpp
//template <typename T> void main_2DKriging(Input in, T context);
void main_2DKriging(Input in, App CoMD);
#endif//CNC||OMP||CIRCLE

#endif //__MAIN_H__
