/** header file for output.cpp
 * **/

#ifndef OUTPUT_HPP
#define OUTPUT_HPP


void printf_calls(int counter, Calls ca);
void printf_timings(int counter, Tms tm);
void printf_colormap_vtk(int i, Node* fields, Input in, int grid_size);
void printf_fields_vtk(int i, Node* node_a, Input in, int grid_size);
void plot_fields(int i, Node* node_a, Input in);

#endif//OUTPUT_HPP
