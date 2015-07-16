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
//  CLASS ImplicitOpenCylinder
//
//=============================================================================



#ifndef ISOEX_IMPLICITOPENCYLINDER_HH
#define ISOEX_IMPLICITOPENCYLINDER_HH

// *** DEBUG ***
#include <iostream>
bool flag_debug_cyl(false);


//== INCLUDES =================================================================

#include "Implicit.hh"
#include <IsoEx/Config/IsoExDefines.hh>

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

	      
/** \class ImplicitOpenCylinder ImplicitOpenCylinder.hh <IsoEx/Implicits/ImplicitOpenCylinder.hh>
    This class implements an implicit open cylinder given by a line and a radius.
    \see IsoEx::Implicit
    \ingroup implicits
*/	      
template< class Vec3 >
class ISOEXDLLEXPORT ImplicitOpenCylinder : public Implicit<Vec3>
{
public:

  typedef typename Vec3::value_type real;
   
  /// \name Constructor & destructor
  //@{

  /// Constructor: given cube center and width
  ImplicitOpenCylinder(const Vec3& _p, const Vec3& _dir, 
                         const real _radius)
    : line_p0_(_p),
      radius_(_radius)
  {
    line_dir_ = _dir;
    line_dir_.normalize();
  }

  /// Empty destructor
  ~ImplicitOpenCylinder() {}

  //@}



  /// \name Abstract interface of implicit objects, see also IsoEx::Implicit.
  //@{

  void compute_dist2line(const Vec3& _point, real & dist) const
  {
    Vec3 v(_point - line_p0_);

    Vec3 v_orth;               // Component of v orthogonal to line_dir_
    v_orth = v - line_dir_ * (v | line_dir_);
    dist = v_orth.norm();
  }

  bool is_inside(const Vec3& _point) const
  {
    real dist;
    compute_dist2line(_point, dist);

    return (dist <= radius_);
  }

  real scalar_distance(const Vec3& _point) const
  {
    real dist;

    compute_dist2line(_point, dist);

    return (dist - radius_);
  }

  bool directed_distance(const Vec3&  _p0,
			 const Vec3&  _p1,
			 Vec3&        _point,
			 Vec3&        _normal,
			 real&        _distance) const
  {
    Vec3 orig(_p0), dirA(_p1-_p0), vA(orig - line_p0_);

    Vec3 dirA_orth;             // Component of dirA orthogonal to line_dir_
    dirA_orth = dirA - line_dir_ * (dirA | line_dir_);
    Vec3 vA_orth;               // Component of vA orthogonal to line_dir_
    vA_orth = vA - line_dir_ * (vA | line_dir_);

    double a = dirA_orth.sqrnorm();
    double b = 2.0*(dirA_orth | vA_orth);
    double c = vA_orth.sqrnorm() - radius_*radius_;
    double d = b*b - 4.0*a*c;

    // *** DEBUG ***
    if (flag_debug_cyl) {
      using namespace std;
      cerr << "*** ImplicitOpenCylinder: radius: " << radius_ << endl;
      cerr << "  Point " << orig;
      cerr << "  dirA: " << dirA << endl;
    }

    if (d >= 0)
    {
      d = sqrt(d);

      double t1 = (-b-d) / (2.0*a);
      double t2 = (-b+d) / (2.0*a);

      // *** DEBUG
      double t;
      if (t1 >= 0.0 && t2 >= 0.0) 
        { t = std::min(t1, t2); }
      else if (t1 > 0.0) 
        { t = t1; }
      else if (t2 > 0.0) 
        { t = t2; }
      else
        { return false ; }

      _point    = orig + dirA*t;

      Vec3 vB = _point - line_p0_;
      Vec3 vB_orth = vB - line_dir_ * (vB | line_dir_);
      _normal   = vB_orth / radius_;
      _distance = ((dirA | _normal) < 0.0) ? dirA.norm()*t : -dirA.norm()*t;

      // *** DEBUG ***
      if (flag_debug_cyl) {
        using namespace std;
        cerr << "  Intersection point: " << _point;
        cerr << "  normal: " << _normal << endl;
        cerr << "  t: " << t
             << "  Distance: " << _distance << endl;
      }

      return true;
    }

    return false;
  }
  
  //@}


protected:

  // Open ylinder is determined by a point, a direction, and a radius.
  // The point and direction define a line.
  // The cylinder is symmetric around the line.
  Vec3  line_dir_;
  Vec3  line_p0_;
  real  radius_;
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_IMPLICITOPENCYLINDER_HH defined
//=============================================================================

