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

//=============================================================================
//
//  CLASS ImplicitCube
//
//=============================================================================



#ifndef ISOEX_IMPLICITCUBE_HH
#define ISOEX_IMPLICITCUBE_HH


//== INCLUDES =================================================================

#include "ImplicitAxisParallelCube.hh"
#include <IsoEx/Config/IsoExDefines.hh>

// *** DEBUG ***
#include <iostream>
bool flag_debug3(false);

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

	      
/** \class ImplicitCube ImplicitCube.hh <IsoEx/Implicits/ImplicitCube.hh>
    This class implements an implicit cube given by its center, width and
    two vectors determining the cube orientation.
    \see IsoEx::Implicit
    \ingroup implicits
*/	      
template< class Vec3 >
class ISOEXDLLEXPORT ImplicitCube : public ImplicitAxisParallelCube<Vec3>
{
public:

  typedef typename Vec3::value_type real;
   
  /// \name Constructor & destructor
  //@{

  /// Constructor: given cube center and width
  ImplicitCube(const Vec3& _center, const Vec3& v0, const Vec3& v1,
               const real _width)
    : ImplicitAxisParallelCube<Vec3>(_center, _width)
  {
    u_[0] = v0;
    u_[0].normalize_cond();
    u_[1] = v1 - u_[0] * (v1 | u_[0]);
    u_[1].normalize_cond();
    u_[2] = u_[0] % u_[1];
  }

  /// Empty destructor
  ~ImplicitCube() {}

  //@}



  /// \name Abstract interface of implicit objects, see also IsoEx::Implicit.
  //@{

  // Get coordinates in frame given by u_[0], u_[1], u_[2].
  void get_cube_frame_coord(const Vec3 & _v_in, Vec3 & _v_out) const
  {
    for (int d = 0; d < 3; d++) 
      { _v_out[d] = (_v_in | u_[d]); }
  }

  // Map coordinates in frame given by u_[0], u_[1], u_[2] to world coord.
  void get_world_coord(const Vec3 & _v_in, Vec3 & _v_out) const
  {
    _v_out = Vec3(0,0,0);
    for (int d = 0; d < 3; d++) 
      { _v_out = _v_out + (_v_in[d] * u_[d]); }
  }

  bool is_inside(const Vec3& _point) const
  {
    Vec3 vA(this->center_ - _point);
    Vec3 vB;

    get_cube_frame_coord(vA, vB);

    return vB.l8_norm() <= this->half_width_;
  }

  real scalar_distance(const Vec3& _point) const
  {
    Vec3 vA(this->center_ - _point);
    Vec3 vB;

    get_cube_frame_coord(vA, vB);

    return vB.l8_norm() - this->half_width_;
  }

  bool directed_distance(const Vec3&  _p0,
			 const Vec3&  _p1,
			 Vec3&        _point,
			 Vec3&        _normal,
			 real&        _distance) const
  {
    Vec3 orig(_p0), dirA(_p1-_p0), vA(orig - this->center_);
    Vec3 dirB, vB;
    Vec3 pointA, normalB;
    real epsilon(0.00001);

    get_cube_frame_coord(dirA, dirB);
    get_cube_frame_coord(vA, vB);

    // *** DEBUG ***
    if (flag_debug3) {
      using namespace std;
      cerr << "*** Point " << orig << endl;
      cerr << "  dirA: " << dirA << "  dirB: " << dirB << endl;
      cerr << "  vA: " << vA << "  vB: " << vB << endl;
    }

    bool is_tmin_set(false), is_tmax_set(false);
    double tmin(0.0), tmax((this->width_+vA.norm())/(dirA.norm()));
    int tmin_inormal(0), tmax_inormal(0);
    for (int d = 0; d < 3; d++) {
      if (dirB[d] > epsilon || dirB[d] < -epsilon) {
        double t0 = (this->half_width_ - vB[d])/dirB[d];
        double t1 = (-this->half_width_ - vB[d])/dirB[d];
        int t0_inormal = 2*d;
        int t1_inormal = 2*d+1;

        if (t0 > t1) {
          std::swap(t0, t1);
          std::swap(t0_inormal, t1_inormal);
        }

        if (t1 < 0) {
          // Ray does not intersect the cube
          return false;
        }
        else {
          if (!is_tmax_set || t1 < tmax) {
            tmax = t1; 
            tmax_inormal = t1_inormal;
            is_tmax_set = true;
          }
        }

        if (t0 > 0) {
          if (!is_tmin_set || t0 > tmin) {
            tmin = t0;
            tmin_inormal = t0_inormal;
            is_tmin_set = true;
          }
        }
      }
      else {
        if (vB[d] > this->half_width_ || vB[d] < -this->half_width_) {
          // ray does not intersect cube
          return false;
        }
      }
    }

    if (is_tmin_set) {
      // Note: If (is_tmin_set == true), then (is_tmax_set == true).

      if (tmin < tmax) {
        _point    = orig + dirA*tmin;
        normalB = this->facet_normal[tmin_inormal];
        get_world_coord(normalB, _normal);
        _distance = ((dirA | _normal) < 0.0) ? 
          dirA.norm()*tmin : -dirA.norm()*tmin;

        // *** DEBUG ***
        if (flag_debug3) {
          using namespace std;
          cerr << "*** Point: " << orig << "  Direction: " << dirA << endl;
          cerr << "  Intersection point: " << _point;
          cerr << "  normal: " << _normal << endl;
          cerr << "  tmin: " << tmin << "  tmax: " << tmax
               << "  Distance: " << _distance << endl;
        }

        return true;
      }
    }
    else if (is_tmax_set) {
      _point    = orig + dirA*tmax;
      normalB = this->facet_normal[tmax_inormal];
      get_world_coord(normalB, _normal);
      _distance = ((dirA | _normal) < 0.0) ? 
        dirA.norm()*tmax : -dirA.norm()*tmax;

      // *** DEBUG ***
      if (flag_debug3) {
        using namespace std;
        cerr << "*** Point: " << orig << "  Direction: " << dirA << endl;
        cerr << "  Intersection point: " << _point;
        cerr << "  normal: " << _normal << endl;
        cerr << "  tmin: " << tmin << "  tmax: " << tmax
             << "  Distance: " << _distance << endl;
      }

      return true;
    }

    return false;
  }
  
  //@}


protected:

  // Cube orientation is determined by three orthogonal unit vectors.
  Vec3  u_[3];

};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_IMPLICITCUBE_HH defined
//=============================================================================

