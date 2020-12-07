/*
  GEO1015.2020
  hw02 
  --
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
*/


#include <iostream>
#include <fstream>
#include <chrono>

//-- to read the params.json file
#include <nlohmann-json/json.hpp>
using json = nlohmann::json;

#include "DatasetASC.h"


//-- struct defined to store the viewpoints
struct Viewpoint3d
{
  double x;
  double y;
  double height;

}; 


//-- forward declarations of functions in this unit, this just ensures that
//-- the compiler knows about the functions
void output_viewshed(DatasetASC &ds, std::vector<Viewpoint3d>& viewpoints, int maxdistance, std::string output_file);


/*
!!! DO NOT MODIFY main() !!!
*/
int main()
{
  // read the params.json JSON file
  std::ifstream i("../data/params.json");
  json j;
  i >> j;
  int maxdistance = j["maxdistance"];
  std::string input_file = j["input_file"];
  std::string output_file = j["output_file"];
  std::vector<Viewpoint3d> viewpoints;
  //-- store each viewpoint in a vector of Viewpoint3d (the struct defined above)
  for (auto each: j["viewpoints"]) {
    Viewpoint3d vp;
    vp.x = each["xy"][0];
    vp.y = each["xy"][1];
    vp.height = each["height"];
    viewpoints.push_back(vp);
  }

  // Record start time
  auto start = std::chrono::high_resolution_clock::now();

  input_file = "../" + input_file; //-- add the parent directory where it is. 
  DatasetASC ds(input_file);
  if (ds.validity == false) {
    return 1;
  }
  // std::cout << ds.ncols << std::endl;
  // std::cout << "The value of the pixel at [row][col]/[0][2] is: " << ds.data[0][2] << std::endl;

  output_viewshed(ds, viewpoints, maxdistance, output_file);

  // Record end time
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << std::fixed << std::setprecision(3) << "--- " << elapsed.count() << " seconds ---" << std::endl;

  //-- we're done, return 0 to say all went fine
  return 0;
}



/*
!!! TO BE COMPLETED !!!
 
Function that calculates the viewshed 
and writes the output raster
 
Input:
    ds            the input dataset (DatasetASC)
    viewpoints:   a list of the viewpoints (x, y, height)
    maxdistance:  max distance one can see
    output_file:  path of the .asc file to write as output
    
Output:
    none (but output .asc file written to 'output-file')
*/
void output_viewshed(DatasetASC &ds, 
                     std::vector<Viewpoint3d>& viewpoints, 
                     int maxdistance, 
                     std::string output_file) {

  std::cout << "==> function output_viewshed() <==" << std::endl;
  std::cout << "      !!! TO BE COMPLETED !!!" << std::endl;

//   int row, col;
//    bool re = ds.xy2rc(505053.0,5258234.5, row, col);
////   bool re = ds.xy2rc(493289,5260558, row, col);
//   std::cout << re << std::endl;
//   std::cout << row << " : " << col << std::endl;
//   double cx, cy;
//   bool re2 = ds.rc2xy(503, 1, cx, cy);
//   std::cout << std::fixed << std::setprecision(1) << "centre" << " : " << cx << ", " << cy << std::endl;


}



