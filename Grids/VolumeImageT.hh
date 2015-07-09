#ifndef _VOLUMEIMAGET_H_
#define _VOLUMEIMAGET_H_

/*******************************************************************************
 * volimage.h
 *
 * Defines a type for a 3D volume image.
 *
 * (C)2005,2006
 * Lehrstuhl I8 RWTH-Aachen, http://www-i8.informatik.rwth-aachen.de
 * Author: Dominik Sibbing
 *
 ******************************************************************************/

//==============================================================================

#include <cassert>

#include <string>
#include <vector>
#include <IsoEx/Grids/ScalarGridT.hh>
#include <IsoEx/Grids/VectorFieldT.hh>
// #include "Types.hh"

// qt
#include <QDateTime>
#include <QDataStream>
#include <QFile>
#include <QFileInfo>
#include <QString>
#include <QDir>

//==============================================================================

#define THRESHOLD 0.05


/** \brief A type for volume images, or 3D textures.
 *
 * The file format is as described in the OpenQVis project. The data consists
 * of a header file and a raw data file. The header file has a format like this:
 *
 * \code
 *  ObjectFileName: CT_Head_large.raw
 *  TaggedFileName: ---
 *  Resolution: 512 512 106
 *  SliceThickness: 0.435547 0.435547 2.0
 *  Origin: 0.0 0.0 0.0
 *  Format: UCHAR
 *  NbrTags: 0
 *  ObjectType: TEXTURE_VOLUME_OBJECT
 *  ObjectModel: RGBA
 *  GridType: EQUIDISTANT
 * \endcode
 *
 * In addition to the OpenQVis format we define two more storage formats,
 * besides UCHAR and USHORT; namely HALF and FLOAT.
 *
 * Also for statistical evaluation, we support an additional fiel called
 * "SigmaData:" which denotes a FLOAT format filename with values that can
 * for example be interpreted as standard variances.
 *
 */

namespace IsoEx
{


template <class Scalar, class Vec3>
class VolumeImageT : public ScalarGridT<Scalar,Vec3>
{
public:


    typedef VectorFieldT<Scalar,3,Vec3> 		VectorField;

    /// Default constructor
    VolumeImageT( const Vec3&  _origin = Vec3( 0,0,0 ),
                  const Vec3&  _x_axis = Vec3( 1,0,0 ),
                  const Vec3&  _y_axis = Vec3( 0,1,0 ),
                  const Vec3&  _z_axis = Vec3( 0,0,1 ),
                  unsigned int            _x_res  = 10,
                  unsigned int            _y_res  = 10,
                  unsigned int            _z_res  = 10 )
            : ScalarGridT<Scalar,Vec3>( _origin, _x_axis, _y_axis, _z_axis, _x_res, _y_res, _z_res )
    {}

    /// Destructor
    virtual ~VolumeImageT() {}


    struct Dimension
    {
        Dimension( unsigned int _x = 1, unsigned int _y = 1, unsigned int _z = 1,
                   Scalar _sx = 1, Scalar _sy = 1, Scalar _sz = 1,
                   Scalar _tx = 0, Scalar _ty = 0, Scalar _tz = 0,
                   Scalar _scale = 1 )
                : x( _x ), y( _y ), z( _z ),
                sx( _sx ), sy( _sy ), sz( _sz ),
                tx( _tx ), ty( _ty ), tz( _tz ),
                scale( _scale )
        { }

        unsigned int   x, y, z;       ///< Dimensions of the image.
        Scalar sx, sy, sz;             ///< Thickness of each dimension.
        Scalar tx, ty, tz;             ///< Translation to origin.
        Scalar scale;                  ///< A scaling factor for the data.

        /** \brief Compares two Dimension structures.
         */
        bool operator!=( const Dimension &di )
        {
            return ( x != di.x &&
                     y != di.y &&
                     z != di.z &&
                     sx != di.sx &&
                     sy != di.sy &&
                     sz != di.sz );
        }

        int max_di(){return std::max( std::max( x,y ),z );}
    };

    struct Neighbor
    {
        Neighbor( int _x, int _y, int _z, Scalar _w ) :
                x( _x ),
                y( _y ),
                z( _z ),
                w_dist( _w )
        {
        }
        int x;
        int y;
        int z;

        Scalar        w_dist;
    };
    typedef std::vector< Neighbor > GaussKernel;


    /// Determines data type of stored image.
    typedef enum Formats
    {
        VI_UCHAR,    ///< Image will use bytes.
        VI_USHORT,   ///< Image will use 16 bit integers.
        VI_HALF,     ///< Image will use 16 bit Scalars.
        VI_FLOAT,    ///< Image will use 32 bit Scalars.
        VI_LOGBYTE   ///< Image will use logarithmic byte encoding.
    } VolumeImageMode;

    void init( const Vec3&  _origin,
               const Vec3&  _x_axis,
               const Vec3&  _y_axis,
               const Vec3&  _z_axis,
               unsigned int            _x_res,
               unsigned int            _y_res,
               unsigned int            _z_res );

    ///load and save
    virtual bool read( const char* _fname );
    virtual bool write( const char* _filename, VolumeImageMode _mode );

    void getByteData( std::vector<unsigned char> &_result ) const;
    void getByteDataGradient( std::vector<unsigned char> &_result ) const;

    void getData( std::vector<float> &_result ) const;
    void getDataGradient( std::vector<float> &_result ) const;

    /// min/max scalar value
    void update_min_max();
    Scalar max_value(){return max_value_;}
    Scalar min_value(){return min_value_;}

    ///Filter
    void bilateral_filter( Scalar _hs, Scalar _hr,  int _iterations );
    void bilateral_filter_simple( Scalar _sigma_s, Scalar _sigma_r,  int _iterations );

    /// normalization and clamping
    void normalize();
    void normalize( double _min, double _max );
    void clamp( double _min, double _max );
    void normalize_with_histogram( Scalar _max, Scalar _percent, int _n_bins = 1000 );
    void histogram( const int _n_bins,
                    std::vector<unsigned int> &_bins,
                    Scalar &_size,
                    bool _to_image = false );

    void set_border();


    ///Gradient calculation
    void assign_gradients();
    void normalize_gradients();
    Vec3 calcGradient( const unsigned int &_x,
                       const unsigned int &_y,
                       const unsigned int &_z );

    /// access to gradient field
    VectorField &grad_field(){return grad_field_;}

private:

    Dimension di_;
    Scalar min_value_;
    Scalar max_value_;

    Scalar threshold_;

    ///Gradients of the VolImage
    VectorField grad_field_;

    ///Functions
    void i_to_xyz( int _i, int &_x, int &_y, int &_z );
    void get_neighbors_in_r( const typename Grid<Vec3>::PointIdx &_p, const int _r, std::vector<unsigned int> &_n );
    void set_bandwidth( Scalar _hs, Scalar _hr );
    void bilateral_update( std::vector<Scalar> &_src,
                           std::vector<Scalar> &_dst,
                           const unsigned int &_point_idx,
                           const Scalar &_hr,
                           const Scalar &_hs );
    void bilateral_update_simple( std::vector<Scalar> &_src,
                                  std::vector<Scalar> &_dst,
                                  int _x, int _y, int _z, Scalar _hs );
    void precalculate_gauss_kernel( Scalar _sigma_s, Scalar _threshold = 0.08 );


    /// Filtering
    GaussKernel cur_gauss_kernel_;
};

} //end namespace IsoEx
#if defined(INCLUDE_TEMPLATES) && !defined(ISOEX_VOLUMEIMAGET_C)
#define ISOEX_VOLUMEIMAGET_TEMPLATES
#include "VolumeImageT.cc"
#endif
// _VOLIMAGE_H_
#endif
