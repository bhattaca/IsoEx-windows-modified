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
#include <IsoEx/Implicits/ImplicitOpenCylinder.hh>
#include <IsoEx/Implicits/ImplicitDoublePlanes.hh>
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
typedef enum { SPHEREII, SPHEREIII, SPHERE_DIFF, CUBE, CUBEII,
               CYLINDER, ANNULUS } 
  LEVEL_SET_TYPE;

//-----------------------------------------------------------------------------


void usage(const char* _argv0)
{
  std::cerr << "\n\nUsage: \n"
            << _argv0 << "  <-e | -m> <-s2 | -s3 | -sdiff | -c | -c2> \n"
            << "      <-a angle> <-r resolution> <-o filename> \n"
            << "      <-center1 x y z> <-center2 x y z>\n"
            << "      <-radius1 r> <-radius2 r>\n"
            << "      <-dir x y z> <-side_dir x y z> \n"
            << "\n";
  
  std::cerr << "  -e      Use Extended Marching Cubes (default)\n"
            << "  -m      Use standard Marching Cubes\n"
            << "  -s2     Union of 2 spheres\n"
            << "  -sdiff  Difference of 2 spheres\n"
            << "  -s3     Combination of 3 spheres (default)\n"
            << "  -c      Cube\n"
            << "  -c2     Union of 2 cubes\n"
            << "  -cyl    Cylinder\n"
            << "  -ann    Annulus\n"
            << "  -a      Feature detection threshold\n"
            << "  -r      Grid resolution (default is 50)\n"
            << "  -o      Write result to filename (should be *.{off,obj,stl}), \n"
            << "          (Default: output.off)\n"
            << "  -center1 x y z  Set center of object 1 to (x,y,z)\n"
            << "  -center1 x y z  Set center of object 1 to (x,y,z)\n"
            << "  -radius1 r      Set radius of sphere 1 to r\n"
            << "  -radius2 r      Set radius of sphere 2 to r\n"
            << "  -dir x y z      Set direction of cube facet to (x,y,z)\n"
            << "  -side_dir x y z Set direction of second cube facet to (x,y,z)\n"
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


int main(int argc, char** argv)
{
  // parameters
  const char*       filename = "output.off";
  unsigned int      res      = 50;
  MODE              mode     = EMC;
  LEVEL_SET_TYPE    object   = SPHEREIII;
  float             angle    = 30.0;
  bool              flag_tilt = true;

  VectorType        center1(0,0,0);
  VectorType        center2(1,0,0);
  VectorType        center3(1,1,0);
  VectorType        dirA(1,0,0);
  VectorType        dirB(0,1,0);
  float             radius1(1.0);
  float             radius2(0.7);
  float             radius3(1.2);
  float             height(1.0);
  float             annulus_width(0.5);
  float             flange_wh(0.2);
  
  
  // parse command line
  int         c;
  extern char *optarg;
  //extern int  optind;
  std::cout <<"argc "<< argc << std::endl;
  int iarg = 1;
  while (iarg < argc) {
   std::cout <<"argv["<< iarg <<"] = "<< argv[iarg] << std::endl;
   std::string s = std::string(argv[iarg]);
   if ( s == "a"){angle = atof(argv[++iarg]);}
   else if ( s == "-e") { mode = EMC;}
   else if ( s == "-m") { mode = MC;}
   else if ( s == "-o") { filename= argv[++iarg];}
   else if ( s == "-r") { res = atoi(argv[++iarg]);}
   else if ( s == "-h") { usage(argv[0]); exit(0);}
   else if ( s == "-s2") { object = SPHEREII; }
   else if ( s == "-s3") { object = SPHEREIII; }
   else if ( s == "-sdiff") { object = SPHERE_DIFF; }
   else if ( s == "-c") { object = CUBE; }
   else if ( s == "-c2") { object = CUBEII; }
   else if ( s == "-no_tilt") { flag_tilt = false; }
   else if ( s == "-cyl") { object = CYLINDER; }
   else if ( s == "-ann") { object = ANNULUS; }
   else if ( s == "-center" || s == "-center1") {
     get_coord(iarg, argc, argv, center1);
     iarg += 3;
   }
   else if ( s == "-center2") {
     get_coord(iarg, argc, argv, center2);
     iarg += 3;
   }
   else if (s == "-radius" || s == "radius1") {
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
		<< "Feature detection angle: " << angle 
		<< std::endl;
      break;
  }
  std::cout << "Grid: " << res << "x" << res << "x" << res << std::endl;
  std::cout << "Output: " << filename << std::endl;

  
  VectorType origin(-2, -2, -2);
  VectorType xaxis(4, 0, 0);
  VectorType yaxis(0, 4, 0);
  VectorType zaxis(0, 0, 4);

  // axis parallel cubes
  ImplicitAxisParallelCube<VectorType>   cubeX1(center1, 1.0);
  ImplicitAxisParallelCube<VectorType>   cubeX2(center2, 1.0);
  CSG::Union<VectorType>     two_cubesX(cubeX1, cubeX2);

  // cubes
  ImplicitCube<VectorType>   cube1(center1, dirA, dirB, 1.0);
  ImplicitCube<VectorType>   cube2(center2, dirA, dirB, 1.0);
  CSG::Union<VectorType>     two_cubes(cube1, cube2);

  // spheres
  ImplicitSphere<VectorType>     sphere1(center1, radius1);
  ImplicitSphere<VectorType>     sphere2(center2, radius2);
  ImplicitSphere<VectorType>     sphere3(center3, radius3);
  CSG::Union<VectorType>         i1(sphere1, sphere2);
  CSG::Difference<VectorType>    i2(i1, sphere3);
  CSG::Difference<VectorType>    diff_s1_s2(sphere1, sphere2);

  // double planes
  ImplicitDoublePlanes<VectorType>   h1(center1, dirA, height);

  // cylinders
  ImplicitOpenCylinder<VectorType> open_cylinder1(center1, dirA, radius1);
  ImplicitOpenCylinder<VectorType> open_cylinder2(center2, dirA, radius2);
  CSG::Intersection<VectorType>    cylinder1(open_cylinder1, h1); 

  // annulus
  ImplicitOpenCylinder<VectorType> 
    open_cylinder3(center1, dirA, radius1+annulus_width/2.0);
  ImplicitOpenCylinder<VectorType> 
    open_cylinder4(center1, dirA, radius1-annulus_width/2.0);
  CSG::Difference<VectorType> 
    diff_ocyl3_ocyl4(open_cylinder3, open_cylinder4);
  CSG::Intersection<VectorType>    annulus1(diff_ocyl3_ocyl4, h1);


  if (object == SPHEREII) {

    ImplicitGrid<VectorType> grid
      (i1, origin, xaxis, yaxis, zaxis, res, res, res);

    extract_mesh(grid, angle, mode, filename);
  }
  else if (object == SPHEREIII) {

    ImplicitGrid<VectorType> grid
      (i2, origin, xaxis, yaxis, zaxis, res, res, res);

    extract_mesh(grid, angle, mode, filename);
  }
  else if (object == SPHERE_DIFF) {
    ImplicitGrid<VectorType> grid
      (diff_s1_s2, origin, xaxis, yaxis, zaxis, res, res, res);

    extract_mesh(grid, angle, mode, filename);
  }
  else if (object == CUBE) {

    if (flag_tilt) {

      ImplicitGrid<VectorType> grid
        (cube1, origin, xaxis, yaxis, zaxis, res, res, res);

      extract_mesh(grid, angle, mode, filename);
    }
    else {
      ImplicitGrid<VectorType> grid
        (cubeX1, origin, xaxis, yaxis, zaxis, res, res, res);

      extract_mesh(grid, angle, mode, filename);
    }
  }
  else if (object == CUBEII) {

    if (flag_tilt) {

      ImplicitGrid<VectorType> grid
        (two_cubes, origin, xaxis, yaxis, zaxis, res, res, res);

      extract_mesh(grid, angle, mode, filename);
    }
    else {
      ImplicitGrid<VectorType> grid
        (two_cubesX, origin, xaxis, yaxis, zaxis, res, res, res);

      extract_mesh(grid, angle, mode, filename);
    }

  }
  else if (object == CYLINDER) {

    ImplicitGrid<VectorType> grid
      (cylinder1, origin, xaxis, yaxis, zaxis, res, res, res);

    extract_mesh(grid, angle, mode, filename);
  }
  else if (object == ANNULUS) {

    ImplicitGrid<VectorType> grid
      (annulus1, origin, xaxis, yaxis, zaxis, res, res, res);

    extract_mesh(grid, angle, mode, filename);
  }
  else {
    std::cerr << "Error. Unknown level set type." << std::endl;
    usage(argv[0]);
    exit(10);
  }

  return 0;
}


//-----------------------------------------------------------------------------
