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
//  CLASS ImplicitAnnulus
//
//=============================================================================



#ifndef ISOEX_IMPLICITANNULUS_HH
#define ISOEX_IMPLICITANNULUS_HH

// *** DEBUG ***
#include <iostream>
bool flag_debug_ann(false);


//== INCLUDES =================================================================

#include "Implicit.hh"
#include <IsoEx/Config/IsoExDefines.hh>

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

	      
/** \class ImplicitAnnulus ImplicitAnnulus.hh <IsoEx/Implicits/ImplicitAnnulus.hh>
    This class implements an implicit closed annulus given by a center, a direction, a radius, a width and a height.
    \see IsoEx::Implicit
    \ingroup implicits
*/	      
template< class Vec3 >
class ISOEXDLLEXPORT ImplicitAnnulus : public Implicit<Vec3>
{
public:

  typedef typename Vec3::value_type real;
   
  /// \name Constructor & destructor
  //@{

  /// Constructor: given annulus center, directionk, radius, width and height
  ImplicitAnnulus(const Vec3& _center, const Vec3& _dir, 
                  const real _radius, const real _width, const real _height)
    : center_(_center),
      radius_(_radius),
      width_(_width),
      height_(_height)
  {
    line_dir_ = _dir;
    line_dir_.normalize();
    half_width_ = _width/2.0;
    half_height_ = _height/2.0;
    inner_radius_ = radius_ - half_width_;
    outer_radius_ = radius_ + half_width_;
  }

  /// Empty destructor
  ~ImplicitAnnulus() {}

  //@}



  /// \name Abstract interface of implicit objects, see also IsoEx::Implicit.
  //@{

  void compute_dist2line(const Vec3& _point, real & dist) const
  {
    Vec3 v(_point - center_);

    Vec3 v_orth;               // Component of v orthogonal to line_dir_
    v_orth = v - line_dir_ * (v | line_dir_);
    dist = v_orth.norm();
  }

  // Compute the (non-negative) distance to the plane through center_
  //   orthogonal to line_dir_
  void compute_dist2plane(const Vec3& _point, real & dist) const
  {
    dist = ((_point - center_) | line_dir_);
    if (dist < 0) { dist = -dist; }
  }

  bool is_inside(const Vec3& _point) const
  {
    real dist2line;
    real dist2plane;
    compute_dist2line(_point, dist2line);
    compute_dist2plane(_point, dist2plane);

    return ( (dist2line >= inner_radius_) &&
             (dist2line <= outer_radius_) && (dist2plane <= half_height_) );
  }

  real scalar_distance(const Vec3& _point) const
  {
    real dist2line;
    real dist2plane;

    compute_dist2line(_point, dist2line);
    compute_dist2plane(_point, dist2plane);

    real dist = (inner_radius_ - dist2line);
    real dist2 = (dist2line - outer_radius_);
    real dist3 = (dist2plane - half_height_);
    if (dist2 > dist) { dist = dist2; }
    if (dist3 > dist) { dist = dist3; }

    // *** DEBUG ***
    if (flag_debug_ann) {
      using namespace std;
      cerr << "In scalar_distance:  Point: " << _point
           << "  dist1: " << inner_radius_-dist2line
           << "  dist2: " << dist2line - outer_radius_
           << "  dist3: " << dist2plane - half_height_
           << endl;
    }

    return (dist);
  }


  // Compute intersection of ray with open inner cylinder.
  // Ray is given by (_orig + _t*_dir).
  bool compute_inner_cylinder_intersection
  (const Vec3&  _orig,
   const Vec3&  _dir,
   real&        _t,
   Vec3&        _point,
   Vec3&        _normal,
   real&        _distance) const
  {
    Vec3 v(_orig - center_);
    real epsilon(0.00001);

    // Initialize
    _t = 0.0;
    _point = Vec3(0.0, 0.0, 0.0);
    _normal = Vec3(0.0, 0.0, 0.0);
    _distance = 0.0;

    // Calculate intersection with annulus
    Vec3 dir_orth;             // Component of _dir orthogonal to line_dir_
    dir_orth = _dir - line_dir_ * (_dir | line_dir_);
    Vec3 v_orth;               // Component of v orthogonal to line_dir_
    v_orth = v - line_dir_ * (v | line_dir_);


    double a = dir_orth.sqrnorm();
    double b = 2.0*(dir_orth | v_orth);
    double c = v_orth.sqrnorm() - inner_radius_*inner_radius_;
    double d = b*b - 4.0*a*c;

    if ((-epsilon <= a) && (a <= epsilon)) { return false; }
    if (d < 0) { return false; }

    d = sqrt(d);

    double t1 = (-b-d) / (2.0*a);
    double t2 = (-b+d) / (2.0*a);

    if (t1 >= 0.0 && t2 >= 0.0) 
      { _t = std::min(t1, t2); }
    else if (t1 >= 0.0) 
      { _t = t1; }
    else if (t2 >= 0.0)
      { _t = t2; }
    else
      { return false; }

    _point    = _orig + _dir*_t;

    Vec3 vB = _point - center_;
    Vec3 vB_orth = vB - line_dir_ * (vB | line_dir_);
    _normal   = -vB_orth / inner_radius_;
    _distance = ((_dir | _normal) < 0.0) ? _dir.norm()*_t : -_dir.norm()*_t;

    return true;
  }

  // Compute intersection of ray with open outer cylinder.
  // Ray is given by (_orig + _t*_dir).
  bool compute_outer_cylinder_intersection
  (const Vec3&  _orig,
   const Vec3&  _dir,
   real&        _t,
   Vec3&        _point,
   Vec3&        _normal,
   real&        _distance) const
  {
    Vec3 v(_orig - center_);
    real epsilon(0.00001);

    // Initialize
    _t = 0.0;
    _point = Vec3(0.0, 0.0, 0.0);
    _normal = Vec3(0.0, 0.0, 0.0);
    _distance = 0.0;

    // Calculate intersection with annulus
    Vec3 dir_orth;             // Component of _dir orthogonal to line_dir_
    dir_orth = _dir - line_dir_ * (_dir | line_dir_);
    Vec3 v_orth;               // Component of v orthogonal to line_dir_
    v_orth = v - line_dir_ * (v | line_dir_);


    double a = dir_orth.sqrnorm();
    double b = 2.0*(dir_orth | v_orth);
    double c = v_orth.sqrnorm() - outer_radius_*outer_radius_;
    double d = b*b - 4.0*a*c;

    if ((-epsilon <= a) && (a <= epsilon)) { return false; }
    if (d < 0) { return false; }

    d = sqrt(d);

    double t1 = (-b-d) / (2.0*a);
    double t2 = (-b+d) / (2.0*a);

    if (t1 >= 0.0 && t2 >= 0.0) 
      { _t = std::min(t1, t2); }
    else if (t1 >= 0.0) 
      { _t = t1; }
    else if (t2 >= 0.0)
      { _t = t2; }
    else
      { return false; }

    _point    = _orig + _dir*_t;

    Vec3 vB = _point - center_;
    Vec3 vB_orth = vB - line_dir_ * (vB | line_dir_);
    _normal   = vB_orth / outer_radius_;
    _distance = ((_dir | _normal) < 0.0) ? _dir.norm()*_t : -_dir.norm()*_t;

    return true;
  }

  // Compute intersection of ray with two planes bounding annulus.
  // Ray is given by (_orig + _t*_dir).
  bool compute_plane_intersection
  (const Vec3&  _orig,
   const Vec3&  _dir,
   real&        _t,
   Vec3&        _point,
   Vec3&        _normal,
   real&        _distance) const
  {
    Vec3 v(_orig - center_);
    real epsilon(0.00001);

    // Initialize
    _t = 0.0;
    _point = Vec3(0.0, 0.0, 0.0);
    _normal = Vec3(0.0, 0.0, 0.0);
    _distance = 0.0;

    real xdir = (_dir | line_dir_);
    real xv = (v | line_dir_);

    if ((-epsilon <= xdir) && (xdir <= epsilon)) { return false; }

    real t1 = (half_height_ - xv) / xdir;
    real t2 = (-half_height_ - xv) / xdir;

    if (xdir > 0) { _normal = line_dir_; }
    else { _normal = -line_dir_; }

    if (t1 > 0.0 && t2 > 0.0) {
      _t = std::min(t1, t2);
      // normal is in opposite direction of _dir.
      _normal = -_normal;
    }
    else if (t1 > 0.0) 
      { _t = t1; }
    else if (t2 > 0.0) 
      { _t = t2; }
    else {
      return false;
    }

    _point    = _orig + _dir*_t;
    _distance = ((_dir | _normal) < 0.0) ? _dir.norm()*_t : -_dir.norm()*_t;

    return true;
  }

  bool directed_distance
  (const Vec3&  _p0,
   const Vec3&  _p1,
   Vec3&        _point,
   Vec3&        _normal,
   real&        _distance) const
  {
    Vec3 orig(_p0), dir(_p1-_p0);
    Vec3 point1, point2, point3, normal1, normal2, normal3;
    real distance1, distance2, distance3;
    real t1, t2, t3;

    // *** DEBUG ***
    if (flag_debug_ann) {
      using namespace std;
      cerr << "*** ImplicitAnnulus: radius: " << radius_ 
           << "  width: " << width_
           << "  height: " << height_
           << "  center: " << center_ << endl;
      cerr << "  Point " << orig;
      cerr << "  dir: " << dir << endl;
    }

    bool ok1 = compute_inner_cylinder_intersection
      (orig, dir, t1, point1, normal1, distance1);
    bool ok2 = compute_outer_cylinder_intersection
      (orig, dir, t2, point2, normal2, distance2);
    bool ok3 = compute_plane_intersection
      (orig, dir, t3, point3, normal3, distance3);


    // *** DEBUG ***
    if (flag_debug_ann) {
      using namespace std;

      real dist2line;
      compute_dist2line(orig, dist2line);

      real dist2plane;
      compute_dist2plane(orig, dist2plane);

      if (ok1) { 
        cerr << "  t1: " << t1 << "  p1: " << point1
             << "  n1 : " << normal1
             << "  d1: " << distance1;
      }
      cerr << "  ok1: " << int(ok1) 
           << endl;
      cerr << "  dist2line: " << dist2line
           << endl;

      if (ok2) { 
        cerr << "  t2: " << t2 << "  p2: " << point2
             << "  n2 : " << normal2
             << "  d2: " << distance2;
      }
      cerr << "  ok2: " << int(ok2) 
           << endl;
      cerr << "  dist2line: " << dist2line
           << endl;


      if (ok3) { 
        cerr << "  t3: " << t3 << "  p3: " << point3
             << "  n3 : " << normal3
             << "  d3: " << distance3;
      }
      cerr << "  ok3: " << int(ok3) 
           << endl;
      cerr << "  dist2plane: " << dist2plane
           << "  Between planes: " << int(dist2plane <= half_height_)
           << endl;
    }


    if (!ok1) {
      real dist2line;
      compute_dist2line(orig, dist2line);
      if (dist2line < inner_radius_) { return false; }
    }

    if (!ok2) {
      real dist2line;
      compute_dist2line(orig, dist2line);
      if (dist2line > outer_radius_) { return false; }
    }

    if (!ok3) {
      real dist2plane;
      compute_dist2plane(orig, dist2plane);
      if (dist2plane > half_height_) { return false; }
    }

    if (ok1) {
      if (ok2) { 
        if (distance1 > distance2) { ok2 = false; }
        else { ok1 = false; }
      }

      if (ok3) { 
        if (distance1 > distance3) { ok3 = false; }
        else { ok1 = false; }
      }
    }

    if (ok2 && ok3) { 
      if (distance2 > distance3) { ok3 = false; } 
      else { ok2 = false; }
    }

    if (ok1) {
      _point = point1;
      _normal = normal1;
      _distance = distance1;

      if (flag_debug_ann)
        { std::cerr << "  Intersect inner cylinder." << std::endl; }

      return true;
    }
    else if (ok2) {
      _point = point2;
      _normal = normal2;
      _distance = distance2;

      if (flag_debug_ann)
        { std::cerr << "  Intersect outer cylinder." << std::endl; }

      return true;
    }
    else if (ok3) {
      _point = point3;
      _normal = normal3;
      _distance = distance3;

      if (flag_debug_ann)
        { std::cerr << "  Intersect bounding planes." << std::endl; }

      return true;
    }

    return false;
  }
  
  //@}


protected:

  // Annulus is determined by a center, a direction, a radius,
  //   a width and a height.
  // The center and direction define a line.
  // The annulus is symmetric around the line.
  Vec3  line_dir_;
  Vec3  center_;
  real  radius_;
  real  inner_radius_;
  real  outer_radius_;
  real  width_;
  real  half_width_;
  real  height_;
  real  half_height_;
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_IMPLICITANNULUS_HH defined
//=============================================================================

