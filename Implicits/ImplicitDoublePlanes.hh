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
//  CLASS ImplicitDoublePlanes
//
//=============================================================================



#ifndef ISOEX_IMPLICITDOUBLEPLANES_HH
#define ISOEX_IMPLICITDOUBLEPLANES_HH


//== INCLUDES =================================================================

#include "Implicit.hh"
#include <IsoEx/Config/IsoExDefines.hh>

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

	      
/** \class ImplicitDoubleplanes ImplicitDoublePlanes.hh <IsoEx/Implicits/ImplicitDoubleplanes.hh>
    This class implements an implicit two planes given by a point, a direction orthogonal to the planes and a distance between the planes.
    \see IsoEx::Implicit
    \ingroup implicits
*/	      
template< class Vec3 >
class ISOEXDLLEXPORT ImplicitDoublePlanes : public Implicit<Vec3>
{
public:

  typedef typename Vec3::value_type real;
   
  /// \name Constructor & destructor
  //@{

  /// Constructor: given cube center and width
  ImplicitDoublePlanes(const Vec3& _p, const Vec3& _orth_dir, 
                    const real _distance_between_planes)
    : center_(_p),
      distance_between_planes_(_distance_between_planes)
  {
    orth_dir_ = _orth_dir;
    orth_dir_.normalize();
    distance_to_center_ = distance_between_planes_ / 2.0;
  }

  /// Empty destructor
  ~ImplicitDoublePlanes() {}

  //@}



  /// \name Abstract interface of implicit objects, see also IsoEx::Implicit.
  //@{

  // Compute the (non-negative) distance to the plane through center_
  void compute_dist2plane(const Vec3& _point, real & dist) const
  {
    dist = ((_point - center_) | orth_dir_);
    if (dist < 0) { dist = -dist; }
  }

  bool is_inside(const Vec3& _point) const
  {
    real dist;
    compute_dist2plane(_point, dist);

    return (dist < distance_to_center_);
  }

  real scalar_distance(const Vec3& _point) const
  {
    real dist;

    compute_dist2plane(_point, dist);

    return (dist - distance_to_center_);
  }

  bool directed_distance(const Vec3&  _p0,
			 const Vec3&  _p1,
			 Vec3&        _point,
			 Vec3&        _normal,
			 real&        _distance) const
  {
    Vec3 orig(_p0), dir(_p1-_p0), v(orig - center_);
    real epsilon(0.00001);


    real xdir = (dir | orth_dir_);
    real xv = (v | orth_dir_);

    if ((-epsilon < xdir) && (xdir < epsilon)) { return false; }

    real t1 = (distance_to_center_ - xv) / xdir;
    real t2 = (-distance_to_center_ - xv) / xdir;

    if (xdir > 0) { _normal = orth_dir_; }
    else { _normal = -orth_dir_; }

    real t;
    if (t1 > 0.0 && t2 > 0.0) {
      t = std::min(t1, t2);

      // normal is in opposite direction of dir.
      _normal = -_normal;
    }
    else if (t1 > 0.0) 
      { t = t1; }
    else if (t2 > 0.0) 
      { t = t2; }
    else {
      return false;
    }

    // *** DEBUG ***
    if (t > 1.00001) { return false; }

    _point    = orig + dir*t;
    _distance = ((dir | _normal) < 0.0) ? dir.norm()*t : -dir.norm()*t;

    return true;
  }
  
  //@}


protected:

  // Double plane is determined by a point, an orthogonal direction, and the distance between planes.
  Vec3  center_;
  Vec3  orth_dir_;
  real  distance_between_planes_;
  real  distance_to_center_;
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_IMPLICITDOUBLEPLANES_HH defined
//=============================================================================

