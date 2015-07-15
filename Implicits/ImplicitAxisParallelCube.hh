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
//  CLASS ImplicitAxisParallelCube
//
//=============================================================================



#ifndef ISOEX_IMPLICITAXISPARALLELCUBE_HH
#define ISOEX_IMPLICITAXISPARALLELCUBE_HH

// *** DEBUG ***
#include <iostream>
bool flag_debug(false);


//== INCLUDES =================================================================

#include "Implicit.hh"
#include <IsoEx/Config/IsoExDefines.hh>

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

	      
/** \class ImplicitAxisParallelCube ImplicitAxisParallelCube.hh <IsoEx/Implicits/ImplicitAxisParellelCube.hh>
    This class implements an implicit axis parallel cube,
    given by its center and width
    \see IsoEx::Implicit
    \ingroup implicits
*/	      
template< class Vec3 >
class ISOEXDLLEXPORT ImplicitAxisParallelCube : public Implicit<Vec3>
{
public:

  typedef typename Vec3::value_type real;
   
  /// \name Constructor & destructor
  //@{

  /// Constructor: given cube center and width
  ImplicitAxisParallelCube(const Vec3& _center, const real _width)
    : center_(_center), 
      width_(_width), 
      half_width_(_width/2.0)
  {
    for (int d = 0; d < 3; d++) {
      facet_normal[2*d] = Vec3(0,0,0);
      facet_normal[2*d][d] = 1;
      facet_normal[2*d+1] = Vec3(0,0,0);
      facet_normal[2*d+1][d] = -1;
    }
  }

  /// Empty destructor
  ~ImplicitAxisParallelCube() {}

  //@}



  /// \name Abstract interface of implicit objects, see also IsoEx::Implicit.
  //@{

  bool is_inside(const Vec3& _point) const
  {
    return (center_ - _point).l8_norm() <= half_width_;
  }

  real scalar_distance(const Vec3& _point) const
  {
    return (center_ - _point).l8_norm() - half_width_;
  }

  bool directed_distance(const Vec3&  _p0,
			 const Vec3&  _p1,
			 Vec3&        _point,
			 Vec3&        _normal,
			 real&        _distance) const
  {
    Vec3 orig(_p0), dir(_p1-_p0), v0(orig - center_);
    real epsilon(0.00001);

    bool is_tmin_set(false), is_tmax_set(false);
    double tmin(0.0), tmax((this->width_+v0.norm())/(dir.norm()));
    int tmin_inormal(0), tmax_inormal(0);
    for (int d = 0; d < 3; d++) {
      if (dir[d] > epsilon || dir[d] < -epsilon) {
        double t0 = (half_width_ - v0[d])/dir[d];
        double t1 = (-half_width_ - v0[d])/dir[d];
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
        if (v0[d] > half_width_ || v0[d] < -half_width_) {
          // ray does not intersect cube
          return false;
        }
      }
    }

    if (is_tmin_set) {
      // Note: If (is_tmin_set == true), then (is_tmax_set == true).

      if (tmin < tmax) {

        _point    = orig + dir*tmin;
        _normal = facet_normal[tmin_inormal];
        _distance = ((dir | _normal) < 0.0) ? 
          dir.norm()*tmin : -dir.norm()*tmin;

        // *** DEBUG ***
        if (flag_debug) {
          using namespace std;
          cerr << "*** Point: " << orig << "  Direction: " << dir << endl;
          cerr << "  Intersection point: " << _point;
          cerr << "  normal: " << _normal << endl;
          cerr << "  tmin: " << tmin << "  tmax: " << tmax
               << "  Distance: " << _distance << endl;
        }

        return true;
      }
    }
    else if (is_tmax_set) {
      _point    = orig + dir*tmax;
      _normal = facet_normal[tmax_inormal];
      _distance = ((dir | _normal) < 0.0) ? 
        dir.norm()*tmax : -dir.norm()*tmax;

      // *** DEBUG ***
      if (flag_debug) {
        using namespace std;
        cerr << "*** Point: " << orig << "  Direction: " << dir << endl;
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

  Vec3  center_;
  real  width_;
  real  half_width_;
  Vec3  facet_normal[6];
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_IMPLICITAXISPARALLELCUBE_HH defined
//=============================================================================

