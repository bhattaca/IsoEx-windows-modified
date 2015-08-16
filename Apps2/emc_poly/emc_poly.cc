// emc_poly.cc
// Created from emc_poly.cc by R. Wenger, 2015.
// Version 0.1.0

/*===========================================================================*\
 *                                                                           *
 *                                IsoEx                                      *
 *        Copyright (C) 2002 by Computer Graphics Group, RWTH Aachen         *
 *                         www.rwth-graphics.de                              *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/


#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <IsoEx/Implicits/ImplicitSphere.hh>
#include <IsoEx/Implicits/ImplicitAxisParallelCube.hh>
#include <IsoEx/Implicits/ImplicitCube.hh>
#include <IsoEx/Implicits/ImplicitCylinder.hh>

#include <IsoEx/Implicits/ImplicitAnnulus.hh>
#include <IsoEx/Implicits/csg.hh>

#include <IsoEx/Grids/ImplicitGridT.hh>

#include <IsoEx/Extractors/MarchingCubesT.hh>
#include <IsoEx/Extractors/ExtendedMarchingCubesT.hh>

//#include <unistd.h>


//#define VectorType OpenMesh::Vec3f
#define VectorType OpenMesh::Vec3d

//-----------------------------------------------------------------------------


using namespace IsoEx;
using namespace OpenMesh;
using namespace OpenMesh::IO;


//-----------------------------------------------------------------------------


// Define the mesh to be used: need vertex normal and status for EMC
struct MyTraits : public DefaultTraits
{
  VertexAttributes   (Attributes::Normal | Attributes::Status);
  HalfedgeAttributes (Attributes::PrevHalfedge);

  /// Use double precision points
  typedef VectorType Point;
  /// Use double precision Normals
  typedef VectorType Normal;
};
typedef TriMesh_ArrayKernelT<MyTraits>  MyMesh;

//-----------------------------------------------------------------------------

typedef enum { MC, EMC } MODE;
typedef enum { SPHERE, CUBE, CYLINDER, ANNULUS } 
  LEVEL_SET_TYPE;
typedef enum { UNION, INTERSECTION, DIFFERENCE } CSG_OPERATION;

//-----------------------------------------------------------------------------


// global variables
VectorType origin(-2, -2, -2);
VectorType xaxis(4, 0, 0);
VectorType yaxis(0, 4, 0);
VectorType zaxis(0, 0, 4);

const char*       filename = "output.off";
unsigned int      res      = 50;
MODE              mode     = EMC;
LEVEL_SET_TYPE    object   = CUBE;
float             feature_angle    = 30.0;


//-----------------------------------------------------------------------------


void usage(const char* _argv0)
{
  std::cerr << "\n\nUsage: \n"
            << _argv0 << "  <-e | -m> <-sphere | -cube | -cyl | -ann> \n"
            << "      <-n N> <-union1> <-union2> <-inter1> <-inter2> <-diff1> <-diff2> \n"
            << "      <-width W> <-height H> <-flange_wh W H>"
            << "      <-center1 x y z> <-center2 x y z>\n"
            << "      <-radius1 r> <-radius2 r>\n"
            << "      <-dir x y z> <-side_dir x y z> \n"
            << "      <-a angle> <-r resolution> <-o filename> <-h>\n"
            << "\n";
  
  std::cerr << "  -e      Use Extended Marching Cubes (default)\n"
            << "  -m      Use standard Marching Cubes\n"
            << "  -sphere Sphere\n"
            << "  -cube   Cube (Default) \n"
            << "  -cyl    Cylinder\n"
            << "  -ann    Annulus\n"
            << "  -n N    Set number of spheres/cubes to N (Default: 1)\n"
            << "  -union1 Set first csg operation to union.\n"
            << "  -union2 Set second csg operation to union.\n"
            << "  -inter1 Set first csg operation to intersection.\n"
            << "  -inter2 Set second csg operation to intersection.\n"
            << "  -diff1  Set first csg operation to difference.\n"
            << "  -diff2  Set second csg operation to difference.\n"
            << "  -width W        Set cube or annulus width to W.\n"
            << "  -height H       Set cylinder or annulus height to H.\n"
            << "  -flange_wh W H  Set flange width to W and height to H\n"
            << "  -center1 x y z  Set center of object 1 to (x,y,z)\n"
            << "  -center2 x y z  Set center of object 2 to (x,y,z)\n"
            << "  -radius1 r      Set radius of sphere 1 or cylinder or annulus to r\n"
            << "  -radius2 r      Set radius of sphere 2 to r\n"
            << "  -dir x y z      Set direction of cube facet to (x,y,z)\n"
            << "  -side_dir x y z Set direction of second cube facet to (x,y,z)\n"
            << "  -a      Feature detection threshold\n"
            << "  -r      Grid resolution (default is 50)\n"
            << "  -o      Write result to filename (should be *.{off,obj,stl}), \n"
            << "          (Default: output.off)\n"
            << "  -version  Print version.\n"
            << "  -h      Print this help message\n"
            << "\n";

  exit(1);
}

//-----------------------------------------------------------------------------

template <typename T>
bool string2val(const char *s, T & x)
{
  std::istringstream x_str;
  
  x_str.str(std::string(s));

  x_str >> x;

  if (!x_str.eof()) 
    { return(false); }
  else
    { return(true); }
}

//-----------------------------------------------------------------------------

void get_coord
(const int iarg, const int argc, char** argv, VectorType & coord)
{
  if (iarg+3 >= argc) {
    std::cerr << "Missing coordinates in argument to " << argv[iarg] 
              << std::endl;
    usage(argv[0]);
    exit(20);
  }

  string2val(argv[iarg+1], coord[0]);
  string2val(argv[iarg+2], coord[1]);
  string2val(argv[iarg+3], coord[2]);
}

//-----------------------------------------------------------------------------


void extract_mesh
(ImplicitGrid<VectorType> & grid, const float angle,
 const  MODE mode, const char * filename)
{
  // extract 0-level isosurface
  MyMesh  mesh;
  switch (mode)
  {
    case MC:  
      grid.build_scalar_distance_cache();
      marching_cubes(grid, mesh); 
      break;

    case EMC:
      grid.build_is_inside_cache();
      extended_marching_cubes(grid, mesh, angle);
      break;
  }

  // write result
  write_mesh(mesh, filename);
}


//-----------------------------------------------------------------------------

// forward declarations
void build_grid_extract_mesh(const Implicit<VectorType> & implicit_);
void extract_from_union(const Implicit<VectorType> & implicit1_,
                        const Implicit<VectorType> & implicit2_);
void extract_from_csg(const Implicit<VectorType> & implicit1_,
                      const Implicit<VectorType> & implicit2_,
                      const CSG_OPERATION csg_op);
void extract_from_csg(const Implicit<VectorType> & implicit1_,
                      const Implicit<VectorType> & implicit2_,
                      const Implicit<VectorType> & implicit3_,
                      const CSG_OPERATION csg_op1,
                      const CSG_OPERATION csg_op2);
void error_too_many_objects();

//-----------------------------------------------------------------------------


int main(int argc, char** argv)
{
  // parameters
  CSG_OPERATION     csg_op1(UNION);
  CSG_OPERATION     csg_op2(DIFFERENCE);
  int               num_objects(1);
  VectorType        center1(0,0,0);
  VectorType        center2(0.5,0.5,0.5);
  VectorType        center3(0.2, 0.2, 1);  
  VectorType        dirA(1,0,0);
  VectorType        dirB(0,1,0);
  VectorType        translation_coord(0.5, 0.5, 0.5);
  float             radius1(1.0);
  float             radius2(0.9);
  float             radius3(0.7);
  float             cube_width(2.0);
  float             height(1.0);
  float             width(1.0);
  float             annulus_width(0.5);
  float             flange_width(0.5);
  float             flange_height(0.5);
  bool              flag_flange(false);
  bool              flag_tilt = true;  
  const std::string VERSION("v0.1.0");
  
  // parse command line
  int         c;
  extern char *optarg;
  //extern int  optind;
  std::cout <<"argc "<< argc << std::endl;
  int iarg = 1;
  while (iarg < argc) {
   std::cout <<"argv["<< iarg <<"] = "<< argv[iarg] << std::endl;
   std::string s = std::string(argv[iarg]);
   if ( s == "a") { feature_angle = atof(argv[++iarg]);}
   else if ( s == "-e") { mode = EMC; }
   else if ( s == "-m") { mode = MC; }
   else if ( s == "-o") { filename= argv[++iarg]; }
   else if ( s == "-r") { res = atoi(argv[++iarg]); }
   else if ( s == "-h") { usage(argv[0]); exit(0); }
   else if ( s == "-sphere") { object = SPHERE; }
   else if ( s == "-cube") { object = CUBE; }
   else if ( s == "-no_tilt") { flag_tilt = false; }
   else if ( s == "-cyl") { object = CYLINDER; }
   else if ( s == "-ann") { object = ANNULUS; }
   else if ( s == "-n") { num_objects = atoi(argv[++iarg]); }
   else if ( s == "-union" || s == "-union1" ) { csg_op1 = UNION; }
   else if ( s == "-inter" || s == "-inter1" ) { csg_op1 = INTERSECTION; }
   else if ( s == "-diff" || s == "-diff1" ) { csg_op1 = DIFFERENCE; }
   else if ( s == "-union2" ) { csg_op2 = UNION; }
   else if ( s == "-inter2" ) { csg_op2 = INTERSECTION; }
   else if ( s == "-diff2" ) { csg_op2 = DIFFERENCE; }
   else if ( s == "-width" ) {
     string2val(argv[++iarg], width);
     cube_width = width;
     annulus_width = width;
   }
   else if ( s == "-height" ) {
     string2val(argv[++iarg], height);
   }
   else if ( s == "-center" || s == "-center1") {
     get_coord(iarg, argc, argv, center1);
     iarg += 3;
   }
   else if ( s == "-center2") {
     get_coord(iarg, argc, argv, center2);
     iarg += 3;
   }
   else if ( s == "-radius" || s == "-radius1") {
     string2val(argv[++iarg], radius1);
   }
   else if (s == "-radius2") {
     string2val(argv[++iarg], radius2);
   }
   else if ( s == "-dir") {
     get_coord(iarg, argc, argv, dirA);
     iarg += 3;
   }
   else if ( s == "-side_dir") {
     get_coord(iarg, argc, argv, dirB);
     iarg += 3;
   }
   else if ( s == "-flange_wh" ) {
     string2val(argv[++iarg], flange_width);
     string2val(argv[++iarg], flange_height);
     flag_flange = true;
   }
   else if ( s == "-version" ) {
     std::cout << "Version: " << VERSION << std::endl;
   }
   else { 
     std::cerr <<"Unrecognised. \n"<< argv[iarg] << std::endl;
     usage( argv[0]);
   }
   iarg++;
  }
  

  // output parameters
  switch (mode)
  {
    case MC:  
      std::cout << "Standard Marching Cubes\n"; 
      break;

    case EMC: 
      std::cout << "Extended Marching Cubes\n"
		<< "Feature detection angle: " << feature_angle 
		<< std::endl;
      break;
  }
  std::cout << "Grid: " << res << "x" << res << "x" << res << std::endl;
  std::cout << "Output: " << filename << std::endl;

  
  if (object == SPHERE) {

    if (num_objects == 1) {
      ImplicitSphere<VectorType>     sphere1(center1, radius1);
      build_grid_extract_mesh(sphere1);
    }
    else if (num_objects == 2) {
      ImplicitSphere<VectorType>     sphere1(center1, radius1);
      ImplicitSphere<VectorType>     sphere2(center2, radius2);

      extract_from_csg(sphere1, sphere2, csg_op1);
    }
    else if (num_objects == 3) {
      ImplicitSphere<VectorType>     sphere1(center1, radius1);
      ImplicitSphere<VectorType>     sphere2(center2, radius2);
      ImplicitSphere<VectorType>     sphere3(center3, radius3);

      extract_from_csg(sphere1, sphere2, sphere3, csg_op1, csg_op2);
    }
    else 
      { error_too_many_objects(); }

  }
  else if (object == CUBE) {

    if (flag_tilt) {

      if (num_objects == 1) {
        ImplicitCube<VectorType>   cube1(center1, dirA, dirB, cube_width);
        build_grid_extract_mesh(cube1);
      }
      else if (num_objects == 2) {
        ImplicitCube<VectorType>   cube1(center1, dirA, dirB, cube_width);
        ImplicitCube<VectorType>   cube2(center2, dirA, dirB, cube_width);;

        extract_from_csg(cube1, cube2, csg_op1);
      }
      else if (num_objects == 3) {
        ImplicitCube<VectorType>   cube1(center1, dirA, dirB, cube_width);;
        ImplicitCube<VectorType>   cube2(center2, dirA, dirB, cube_width);;
        ImplicitCube<VectorType>   cube3(center3, dirA, dirB, cube_width);;

        extract_from_csg(cube1, cube2, cube3, csg_op1, csg_op2);
      }
      else 
        { error_too_many_objects(); }
    }
    else {

      if (num_objects == 1) {
        ImplicitAxisParallelCube<VectorType>   cubeX1(center1, cube_width);;
        build_grid_extract_mesh(cubeX1);
      }
      else if (num_objects == 2) {
        ImplicitAxisParallelCube<VectorType>   cubeX1(center1, cube_width);;
        ImplicitAxisParallelCube<VectorType>   cubeX2(center2, cube_width);;

        extract_from_csg(cubeX1, cubeX2, csg_op1);
      }
      else if (num_objects == 3) {
        ImplicitAxisParallelCube<VectorType>   cubeX1(center1, cube_width);;
        ImplicitAxisParallelCube<VectorType>   cubeX2(center2, cube_width);;
        ImplicitAxisParallelCube<VectorType>   cubeX3(center3, cube_width);;

        extract_from_csg(cubeX1, cubeX2, cubeX3, csg_op1, csg_op2);
      }
      else 
        { error_too_many_objects(); }
    }
  }
  else if (object == CYLINDER) {

    if (flag_flange) {
      ImplicitCylinder<VectorType> 
        cylinder1(center1, dirA, radius1+2*flange_width, height);
      ImplicitCylinder<VectorType> 
        cylinder2(center1, dirA, radius1, height+2*flange_height);

      extract_from_union(cylinder1, cylinder2);
    }
    else {
      ImplicitCylinder<VectorType> cylinder1(center1, dirA, radius1, height);
      build_grid_extract_mesh(cylinder1);
    }

  }
  else if (object == ANNULUS) {

    if (flag_flange) {
      ImplicitAnnulus<VectorType> 
        annulus1(center1, dirA, radius1, annulus_width+2*flange_width, height);
      ImplicitAnnulus<VectorType> 
        annulus2(center1, dirA, radius1, annulus_width, height+2*flange_height);
      
      extract_from_union(annulus1, annulus2);
    }
    else {
      ImplicitAnnulus<VectorType> 
        annulus1(center1, dirA, radius1, annulus_width, height);
      build_grid_extract_mesh(annulus1);
    }

  }
  else {
    std::cerr << "Error. Unknown level set type." << std::endl;
    usage(argv[0]);
    exit(10);
  }

  return 0;
}


//-----------------------------------------------------------------------------

// Build implicit grid and then extract mesh.
void build_grid_extract_mesh(const Implicit<VectorType> & implicit_)
{
  ImplicitGrid<VectorType> 
    grid(implicit_, origin, xaxis, yaxis, zaxis, res, res, res);

  extract_mesh(grid, feature_angle, mode, filename);
}

void extract_from_union(const Implicit<VectorType> & implicit1_,
                        const Implicit<VectorType> & implicit2_)
{
  CSG::Union<VectorType>     union_1_2(implicit1_, implicit2_);
  build_grid_extract_mesh(union_1_2);
}


void extract_from_intersection(const Implicit<VectorType> & implicit1_,
                               const Implicit<VectorType> & implicit2_)
{
  CSG::Intersection<VectorType>   intersection_1_2(implicit1_, implicit2_);
  build_grid_extract_mesh(intersection_1_2);
}


void extract_from_difference(const Implicit<VectorType> & implicit1_,
                             const Implicit<VectorType> & implicit2_)
{
  CSG::Difference<VectorType>   difference_1_2(implicit1_, implicit2_);
  build_grid_extract_mesh(difference_1_2);
}

void undefined_csg_operation()
{
  using namespace std;
  cerr << "Programming error.  Undefined CSG operation." << endl;
  cerr << "Exiting..." << endl;
  exit(200);
}

void extract_from_csg(const Implicit<VectorType> & implicit1_,
                      const Implicit<VectorType> & implicit2_,
                      const CSG_OPERATION csg_op)
{
  if (csg_op == UNION) 
    { extract_from_union(implicit1_, implicit2_); }
  else if (csg_op == INTERSECTION) 
    { extract_from_intersection(implicit1_, implicit2_); }
  else if (csg_op == DIFFERENCE) 
    { extract_from_difference(implicit1_, implicit2_); }
  else
    { undefined_csg_operation(); }
}

void extract_from_csg(const Implicit<VectorType> & implicit1_,
                      const Implicit<VectorType> & implicit2_,
                      const Implicit<VectorType> & implicit3_,
                      const CSG_OPERATION csg_op1,
                      const CSG_OPERATION csg_op2)
{
  if (csg_op1 == UNION) {
    CSG::Union<VectorType>     union_1_2(implicit1_, implicit2_);
    extract_from_csg(union_1_2, implicit3_, csg_op2);
  }
  else if (csg_op1 == INTERSECTION) {
    CSG::Intersection<VectorType>   intersection_1_2(implicit1_, implicit2_);
    extract_from_csg(intersection_1_2, implicit3_, csg_op2);
  }
  else if (csg_op1 == DIFFERENCE) {
    CSG::Difference<VectorType>   difference_1_2(implicit1_, implicit2_);
    extract_from_csg(difference_1_2, implicit3_, csg_op2);
  }
  else
    { undefined_csg_operation(); }
}

void error_too_many_objects()
{
  using namespace std;
  cerr << "Too many objects.  Reduce argument of -n." << endl;
  cerr << "Exiting..." << endl;
  exit(210);
}
