/*******************************************************************************
 * volimage.cc
 *
 * Defines a type for a 3D volume image.
 *
 * (C)2005,2006
 * Lehrstuhl I8 RWTH-Aachen, http://www-i8.informatik.rwth-aachen.de
 * Author: Dominik Sibbing
 *
 ******************************************************************************/

#define ISOEX_VOLUMEIMAGET_C

//==============================================================================

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include <limits>

#include <QImage>

#include "VolumeImageT.hh"

namespace IsoEx
{


template <class Scalar, class Vec3>
bool
VolumeImageT<Scalar,Vec3>::read( const char* _fname )
{
    // get file infos
    QFileInfo finfo( _fname );

    QDir    dir( finfo.absolutePath() );
    QString base_name = finfo.baseName();
    QString suffix = finfo.suffix();

    QString header_name = finfo.absoluteFilePath();
    QString rawfile_name = "";

    // mode and dimension
    VolumeImageMode mode = VI_UCHAR;
    Dimension di;

    if ( suffix == "hdr" )
    {
        //      - .hdr File �ffnen
        // - int32_t einlesen -> HDR length
        //  - Falls HDR length != 384 -> abort
        // - Seek an Position 0x28 vom Dateianfang an
        // - int16_t dim[8] einlesen
        // - dim[0] == 3 -> anzahl dimensionen
        // - dim[1] bis dim[3] sind X Y Z
        // - seek an position 0x46 vom Dateianfang an
        // - int16_t datatype einlesen
        // - datatype == 4 ist short (rest in .h file erkl�rt)
        //   -> deren short ist signed -- vorsicht, bzw. volimage erweitern!
        // - volimage datentyp setzen
        // - volimage puffer auf X Y Z setzen
        // - .img file �ffnen und einlesen

        QFile header( header_name );
        if ( header.open( QIODevice::ReadOnly ) )
        {
            QDataStream stream( &header );
            stream.setByteOrder( QDataStream::LittleEndian );

            // check length
            qint32 buf32;
            stream >> buf32;
            if ( buf32 != 348 )
            {
                std::cout << " ERROR: wrong file length: "
                        << header_name.toUtf8().data() << ": " << buf32
                        << std::endl;
                return false;
            }

            // seek to position
            stream.skipRawData( 0x24 );

            // read dimension
            qint16 buf16;
            stream >> buf16;
            if ( buf16 > 4 )
            {
                std::cout << " ERROR: dimension != 3. " << ": " << buf16 << std::endl;
                return false;
            }

            stream >> buf16; di.x = buf16;
            stream >> buf16; di.y = buf16;
            stream >> buf16; di.z = buf16;

            // data type
            stream.skipRawData( 22 );
            stream >> buf16;

            switch ( buf16 )
            {
                case 2:
                    mode = VI_UCHAR;
                    break;
                case 4:
                    mode = VI_USHORT;
                    break;
                case 16:
                    mode = VI_FLOAT;
                    break;
                default:
                    std::cout << " ERROR: Unsopported data format." << std::endl;
                    return false;
            }

            // set rawfile
            QFileInfo tmp( dir, base_name + ".img" );
            rawfile_name = tmp.absoluteFilePath();
        }
        else
        {
            std::cout << "ERROR: reading header file: "
                    << header_name.toUtf8().data() << std::endl;
            return false;
        }
    }
    else if ( suffix == "dat" )
    {
        std::ifstream header( header_name.toUtf8().data(), std::ios::in );
        if ( !header.is_open() )
        {
            std::cout << "ERROR: reading header file." << std::endl;
            return false;
        }

        std::string buf;

        while ( !header.eof() )
        {
            header >> buf;
            if ( buf == "ObjectFileName:" )
            {
                header >> buf;

                QFileInfo tmp( dir, QString( buf.data() ) );
                rawfile_name = tmp.absoluteFilePath();
            }
            else if ( buf == "Resolution:" )
            {
                header >> di.x >> di.y >> di.z;
            }
            else if ( buf == "SliceThickness:" )
            {
                header >> di.sx >> di.sy >> di.sz;
            }
            else if ( buf == "Origin:" )
            {
                header >> di.tx >> di.ty >> di.tz;
            }
            else if ( buf == "Format:" )
            {
                header >> buf;
                if ( buf == "UCHAR" )
                    mode = VI_UCHAR;
                else if ( buf == "USHORT" )
                    mode = VI_USHORT;
                else if ( buf == "HALF" )
                    mode = VI_HALF;
                else if ( buf == "FLOAT" )
                    mode = VI_FLOAT;
                else if ( buf == "LOGBYTE" )
                    mode = VI_LOGBYTE;
            }
            else if ( buf == "Scale:" )
            {
                header >> di.scale;
            }
        }

        header.close();
    }
    else
    {
        std::cout << "ERROR: unrecognizable file format. " << std::endl;
        return false;
    }


    Vec3 origin( di.tx, di.ty, di.tz );
    Vec3 x_axis( Scalar( di.x )*di.sx,0.0,0.0 );
    Vec3 y_axis( 0.0,Scalar( di.y )*di.sy,0.0 );
    Vec3 z_axis( 0.0,0.0,Scalar( di.z )*di.sz );

    this->initialize( origin, x_axis, y_axis, z_axis, di.x, di.y, di.z );
    grad_field_.initialize( origin, x_axis, y_axis, z_axis, di.x, di.y, di.z );

    di_ = di;

    std::ifstream datafile( rawfile_name.toUtf8().data(), std::ios::binary );
    if ( !datafile.is_open() )
    {
        std::cout << "ERROR: reading data file: "
                << rawfile_name.toUtf8().data() << std::endl;
        return false;
    }
    int x = 0;
    int y = 0;
    int z = 0;

    min_value_ =  std::numeric_limits<Scalar>::max();
    max_value_ = -std::numeric_limits<Scalar>::max();


    switch ( mode )
    {

        case VI_UCHAR:
            for ( unsigned int i = 0; i < di.x * di.y * di.z; ++i )
            {
                unsigned char item;
                datafile.read( reinterpret_cast<char *>( &item ) , sizeof( unsigned char ) );

                i_to_xyz( i,x,y,z );
                this->value( x,y,z ) = di.scale * item / 255.0f;

                if ( this->value( x,y,z ) > max_value_ ) max_value_ = this->value( x,y,z );
                if ( this->value( x,y,z ) < min_value_ ) min_value_ = this->value( x,y,z );
            }
            break;

        case VI_USHORT:
            for ( unsigned int i = 0; i < di.x * di.y * di.z; ++i )
            {
                unsigned short item;
                datafile.read( reinterpret_cast<char *>( &item ), sizeof( unsigned short ) );
                i_to_xyz( i,x,y,z );
                this->value( x,y,z ) = di.scale * item / 65535.0f;

                if ( this->value( x,y,z ) > max_value_ ) max_value_ = this->value( x,y,z );
                if ( this->value( x,y,z ) < min_value_ ) min_value_ = this->value( x,y,z );
            }
            break;

        case VI_FLOAT:
            for ( unsigned int i = 0; i < di.x * di.y * di.z; ++i )
            {
                float item;
                datafile.read( reinterpret_cast<char *>( &item ), sizeof( float ) );
                i_to_xyz( i,x,y,z );
                this->value( x,y,z ) = di.scale *Scalar( item );

                if ( this->value( x,y,z ) > max_value_ ) max_value_ = this->value( x,y,z );
                if ( this->value( x,y,z ) < min_value_ ) min_value_ = this->value( x,y,z );
            }
            break;

        case VI_HALF:
            for ( unsigned int i = 0; i < di.x * di.y * di.z; ++i )
            {
                Scalar item;
                datafile.read( reinterpret_cast<char *>( &item ), sizeof( Scalar ) );
                i_to_xyz( i,x,y,z );
                this->value( x,y,z ) = di.scale * item;

                if ( this->value( x,y,z ) > max_value_ ) max_value_ = this->value( x,y,z );
                if ( this->value( x,y,z ) < min_value_ ) min_value_ = this->value( x,y,z );
            }
            break;

        case VI_LOGBYTE:
            for ( unsigned int i = 0; i < di.x * di.y * di.z; ++i )
            {
                unsigned char item;
                datafile.read( reinterpret_cast<char *>( &item ),
                               sizeof( unsigned char ) );

                if ( item != 0 )
                {
                    i_to_xyz( i,x,y,z );
                    this->value( x,y,z ) = di.scale * powf( 10, ( 255 - item ) / -25.0f );
                }
                else
                {
                    i_to_xyz( i,x,y,z );
                    this->value( x,y,z ) = 0.0f;
                }

                if ( this->value( x,y,z ) > max_value_ ) max_value_ = this->value( x,y,z );
                if ( this->value( x,y,z ) < min_value_ ) min_value_ = this->value( x,y,z );
            }
            break;

        default:
            std::cout << "VolumeImage::load: Unknown format." << std::endl;
            return false;
            break;

    }

    datafile.close();

    std::cout << "Volume loaded: " << std::endl;
    std::cout << "\tresolution: " <<  this->x_resolution()
    << ", " << this->y_resolution()
    << ", " << this->z_resolution() << std::endl;
    std::cout << "\tvoxelsize : " <<  this->dx().norm()
    << ", " << this->dy().norm()
    << ", " << this->dz().norm() << std::endl;

    return true;
}

//-----------------------------------------------------------------------------

template <class Scalar, class Vec3>
bool
VolumeImageT<Scalar,Vec3>::write( const char* _fname, VolumeImageMode _mode )
{

//     get file infos
//     QFileInfo finfo( _fname );

//     QString fname = finfo.absoluteFilePath() + "/" + finfo.baseName() + ".dat";
//     QString fdata = finfo.absoluteFilePath() + "/" + finfo.baseName() + ".raw";

    std::string fbasename( _fname );

    std::string pathname;
    if ( fbasename[0] == '/' || ( fbasename[0] == '.' && fbasename[1] == '.' ) )
    {
        // filename begins with a path description
        pathname = fbasename.substr( 0, fbasename.find_last_of( '/' ) + 1 );
    }

    std::string filename;
    if ( fbasename[0] == '/' || ( fbasename[0] == '.' && fbasename[1] == '.' ) )
    {
        // filename begins with a path description
        filename = fbasename.substr( fbasename.find_last_of( '/' ) + 1, fbasename.length() - fbasename.find_last_of( '/' ) - 1 );
    }

    const std::string fname   = fbasename + ".dat";
    const std::string fdata   = fbasename + ".raw";

    std::ofstream header( fname.data(), std::ios::out );
    if ( !header.is_open() )
        return false;

    header << "ObjectFileName: " << filename + ".raw" << std::endl;
    header << "TaggedFileName: --- " << std::endl;
    header << "Resolution: " << di_.x << " " << di_.y << " " << di_.z << std::endl;
    header << "SliceThickness: " << di_.sx << " " << di_.sy << " " << di_.sz << std::endl;
    header << "Origin: " << di_.tx << " " << di_.ty << " " << di_.tz << std::endl;
    header << "Format: ";
    switch ( _mode )
    {

        case VI_UCHAR:
            header << "UCHAR";
            break;

        case VI_USHORT:
            header << "USHORT";
            break;

        case VI_FLOAT:
            header << "FLOAT";
            break;

        case VI_HALF:
            header << "HALF";
            break;


        case VI_LOGBYTE:
            header << "LOGBYTE";
            break;

        default:
            std::cout << "VolumeImage::save: Unknown format." << std::endl;
            return false;
            break;

    }
    header << std::endl;
    header << "Scale: " << di_.scale << std::endl;
    header << "NbrTags: 0" << std::endl;
    header << "ObjectType: TEXTURE_VOLUME_OBJECT" << std::endl;
    header << "ObjectModel: RGBA" << std::endl;
    header << "GridType: EQUIDISTANT" << std::endl;

    header.close();

    std::ofstream datafile( fdata.data(), std::ios::binary );
    if ( !datafile.is_open() )
        return false;

    int x = 0;
    int y = 0;
    int z = 0;
    switch ( _mode )
    {

        case VI_UCHAR:
            for ( unsigned int i = 0; i < di_.x * di_.y * di_.z; ++i )
            {
                i_to_xyz( i,x,y,z );

                unsigned char item;
                if ( this->value( x,y,z ) <= 1.0f )
                {
                    item = static_cast<unsigned char>( this->value( x,y,z ) * 255 );
                }
                else
                {
                    item = 255;
                }
                //         if (item > 0) std::cout << (int)item << std::endl;
                datafile.write( reinterpret_cast<char *>( &item ),
                                sizeof( unsigned char ) );
            }
            break;

        case VI_USHORT:
            for ( unsigned int i = 0; i < di_.x * di_.y * di_.z; ++i )
            {
                i_to_xyz( i,x,y,z );
                // HACK u_int16_t item;
                unsigned short item;
                if ( this->value( x,y,z ) <= 1.0f )
                {
                    item = static_cast<unsigned char>( this->value( x,y,z ) * 65535 );
                }
                else
                {
                    item = 65535;
                }
                datafile.write( reinterpret_cast<char *>( &item ),
                                sizeof( unsigned short ) );
            }
            break;

        case VI_FLOAT:
            for ( unsigned int i = 0; i < di_.x * di_.y * di_.z; ++i )
            {
                i_to_xyz( i,x,y,z );
                float item( this->value( x,y,z ) );
                datafile.write( reinterpret_cast<char *>( &item ),
                                sizeof( float ) );
            }
            break;

        case VI_HALF:
            for ( unsigned int i = 0; i < di_.x * di_.y * di_.z; ++i )
            {
                i_to_xyz( i,x,y,z );
                Scalar item( this->value( x,y,z ) );
                datafile.write( reinterpret_cast<char *>( &item ),
                                sizeof( Scalar ) );
            }
            break;

        case VI_LOGBYTE:
        {
            const Scalar cf = powf( 10.0f, -10.24f );

            for ( unsigned int i = 0; i < di_.x * di_.y * di_.z; ++i )
            {
                i_to_xyz( i,x,y,z );
                unsigned char item;

                if ( this->value( x,y,z ) <= 1.0f && this->value( x,y,z ) >= cf )
                {
                    item = static_cast<unsigned char>( 255 - log10f( this->value( x,y,z ) ) * -25 );
                }
                else if ( this->value( x,y,z ) >= 1.0f )
                {
                    item = 255;
                }
                else
                {
                    item = 0;
                }

                datafile.write( reinterpret_cast<char *>( &item ),
                                sizeof( unsigned char ) );
            }
        }
        break;

    }

    datafile.close();

    std::cout << "Image wrote to " << _fname << std::endl;
    return true;
}

//-----------------------------------------------------------------------------

template <class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::update_min_max()
{
    unsigned int x_res( this->x_resolution() );
    unsigned int y_res( this->y_resolution() );
    unsigned int z_res( this->z_resolution() );

    min_value_ =  std::numeric_limits<Scalar>::max();
    max_value_ = -std::numeric_limits<Scalar>::max();

    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                if ( this->value( x,y,z ) > max_value_ ) max_value_ = this->value( x,y,z );
                if ( this->value( x,y,z ) < min_value_ ) min_value_ = this->value( x,y,z );
            }

    std::cout << "Min/Max value: " << min_value_ << " / " << max_value_ << std::endl;
}


//-----------------------------------------------------------------------------


/** \brief Returns an array of bytes.
*
* Returns the volume data converted into bytes, e.g. for usage as a
* 3D texture. The range 0.0 - 1.0 is converted to 0 - 255. Values that are
* out of these bounds are clipped to the nearest bound.
*/
template <class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::getByteData( std::vector<unsigned char> &_result ) const
{

    unsigned int x_res( this->x_resolution() );
    unsigned int y_res( this->y_resolution() );
    unsigned int z_res( this->z_resolution() );

    _result.clear();
    _result.resize( x_res*y_res*z_res );
    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                unsigned int i = x_res*y_res*z + x_res*y + x;
                if ( this->value( x,y,z ) < 0.0f )
                    _result[i] = 0;
                else if ( this->value( x,y,z ) > 1.0f )
                    _result[i] = 255;
                else
                    _result[i] = ( unsigned char )( this->value( x,y,z ) * 255 );
            }
    std::cout << "Byte data ready" << std::endl;
}


//-----------------------------------------------------------------------------


template <class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::getData( std::vector<float> &_result ) const
{
    unsigned int x_res( this->x_resolution() );
    unsigned int y_res( this->y_resolution() );
    unsigned int z_res( this->z_resolution() );

    _result.clear();
    _result.resize( x_res*y_res*z_res );
    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                unsigned int i = x_res*y_res*z + x_res*y + x;
                _result[i] = this->value( x,y,z );
            }
    std::cout << "Data ready" << std::endl;
}


//-----------------------------------------------------------------------------


template <class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::i_to_xyz( int _i, int &_x, int &_y, int &_z )
{
    _x = _i % di_.x; _i /= di_.x;
    _y = _i % di_.y; _i /= di_.y;
    _z = _i;
}


//-----------------------------------------------------------------------------


template <class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::get_neighbors_in_r( const typename Grid<Vec3>::PointIdx &_p,
        const int _r,
        std::vector<unsigned int> &_n )
{
    static int x;
    static int y;
    static int z;
    i_to_xyz( _p, x,y,z );

    static int nx;
    static int ny;
    static int nz;

    static int a;
    static int b;

    static int j;
    static int i;
    static int k;

    _n.clear();
    for ( j = -_r; j <= _r; ++j )
    {
        a = int( sqrt( Scalar( _r*_r ) - Scalar( j*j ) ) );
        for ( i = -a; i <= a; ++i )
        {
            b = int( sqrt( Scalar( a*a ) - Scalar( i*i ) ) );
            for ( k = -b; k <= b; ++k )
            {
                nx = x+k;
                ny = y+i;
                nz = z+j;
                if ( nx >= 0 && nx < ( int )di_.x && ny >= 0 && ny < ( int )di_.y && nz >= 0 && nz < ( int )di_.z )
                {
                    _n.push_back( nz*di_.x*di_.y+ny*di_.x+nx );
                    //             if (pointNode)
                    //             {
                    //               pointNode->add_color(ACG::Vec3uc(255,0,0));
                    //               pointNode->add_point(this->point(nz*di_.x*di_.y+ny*di_.x+nx));
                    //             }
                }
            }
        }
    }
}



//   ///Bilateral Filtering
//   //=============================================================================
//

template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::
precalculate_gauss_kernel( Scalar _sigma_s, Scalar _threshold )
{
    // calculate radius of neighbors ( dist > 3*sigma -> weight < 0.02 )
    int r = int( _sigma_s * 4.0 + 0.5 );

    Scalar pre_w = 1.0/( _sigma_s * sqrt( 2.0 * M_PI ) );

    cur_gauss_kernel_.clear();
    for ( int z = -r; z <= r; ++z )
        for ( int y = -r; y <= r; ++y )
            for ( int x = -r; x <= r; ++x )
            {
                // distance to local origin
//             Vec3 p = this->point( x, y, z );
//             Scalar dist2 = ( ( p - this->origin() ) | ( p - this->origin() ) );
                Scalar dist = x*x + y*y + z*z;

                // weight
                Scalar w = pre_w * exp( -( dist / ( 2.0*_sigma_s*_sigma_s ) ) );

                // set new entry
                if ( w > _threshold )
                    cur_gauss_kernel_.push_back( Neighbor( x, y, z, w ) );
            }

    std::cout << cur_gauss_kernel_.size() << " Neighbors added to Kernel." << std::endl;
}

//=============================================================================

template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::bilateral_filter( Scalar _hs, Scalar _hr,  int _iterations )
{
    QTime stop_watch;
    stop_watch.start();

    //precalculation of the ranges
    set_bandwidth( _hs, _hr );

    std::vector<Scalar> tmp1;
    std::vector<Scalar> tmp2;
    tmp1.clear();
    tmp2.clear();

    //copy data of volimage into x
    for ( unsigned int z = 0; z < di_.z; ++z )
        for ( unsigned int y = 0; y < di_.y; ++y )
            for ( unsigned int x = 0; x < di_.x; ++x )
                tmp1.push_back( this->value( x,y,z ) );

    tmp2.resize( tmp1.size(),0.0 );

    int i = 0;
    int mod = this->n_points() / 1000;
    for ( ; i < _iterations; ++i )
    {
// #pragma omp parallel for
        for ( unsigned int j = 0; j < this->n_points() ; ++j )
        {
            if ( i%2==0 ) bilateral_update( tmp1,tmp2,j,_hr,_hs );
            else        bilateral_update( tmp2,tmp1,j,_hr,_hs );

            if ( j % mod == 0 )
                std::cout  << j << " | " << this->n_points() << "                           \r"; std::cout.flush();

        }
        std::cout << "Iteration: " << i << std::endl;
    }

    //store new pixeldata
    for ( unsigned int z = 0; z < di_.z; ++z )
        for ( unsigned int y = 0; y < di_.y; ++y )
            for ( unsigned int x = 0; x < di_.x; ++x )
            {
                if ( i%2==0 )
                {
                    this->value( x,y,z ) = tmp1[di_.x*di_.y*z + di_.x*y + x];
                }
                else
                {
                    this->value( x,y,z ) = tmp2[di_.x*di_.y*z + di_.x*y + x];
                }
            }
    //reset ranges

    set_bandwidth( 1.0/_hs, 1.0/_hr );

    unsigned int elapsed = stop_watch.elapsed();
    std::cout << "Time for filtering: " << elapsed << "ms" << std::endl;
}

//=============================================================================

template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::bilateral_update( std::vector<Scalar> &_src,
                                        std::vector<Scalar> &_dst,
                                        const unsigned int &_point_idx,
                                        const Scalar &_hr,
                                        const Scalar &_hs )
{
    std::vector<unsigned int> neighbors;
    neighbors.clear();

    get_neighbors_in_r( _point_idx, int( _hs )+1, neighbors );

    //     std::cout << neighbors.size() << std::endl;

    Scalar acc_w = 0.0;
    Scalar I = 0.0;

    Vec3 pos( this->point( _point_idx ) );

    Vec3 p_pn;

    if ( pos[0] >= Scalar( di_.x-3 )/_hs || pos[0] <= 2.0/_hs ||
            pos[1] >= Scalar( di_.y-3 )/_hs || pos[1] <= 2.0/_hs ||
            pos[2] >= Scalar( di_.z-3 )/_hs || pos[2] <= 2.0/_hs ||
            _src[_point_idx] < threshold_ )
    {
        _dst[_point_idx] = 0.0;
        //       if (_src[_point_idx] < threshold_) std::cout << "Test!" << std::endl;
        return;
    }


    Scalar dist = 0.0;
    Scalar weight = 0.0;

    //     int i = 0;

    for ( unsigned int iter = 0; iter < neighbors.size(); ++iter )
    {
        p_pn = pos - this->point( neighbors[iter] );

        dist = ( p_pn|p_pn ) + ( _src[_point_idx]-_src[neighbors[iter]] )*( _src[_point_idx]-_src[neighbors[iter]] );
        //       std::cout << dist << std::endl;

        weight = exp( -0.5*dist );

        acc_w  += weight;
        I += _src[neighbors[iter]]*weight;
    }
    //     std::cout << i << ", " << neighbors.size()<<std::endl;

    I /= acc_w;
    //     if (I > 0.1) std::cout << I << std::endl;
    _dst[_point_idx] = I;
}

//=============================================================================

template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::bilateral_filter_simple( Scalar _sigma_s, Scalar _sigma_r,  int _iterations )
{
    QTime stop_watch;
    stop_watch.start();

    // precalculate gauss kernel
    precalculate_gauss_kernel( _sigma_s );

    std::vector<Scalar> tmp1;
    std::vector<Scalar> tmp2;
    tmp1.clear();
    tmp2.clear();

    //copy data of volimage into x
    for ( unsigned int z = 0; z < di_.z; ++z )
        for ( unsigned int y = 0; y < di_.y; ++y )
            for ( unsigned int x = 0; x < di_.x; ++x )
                tmp1.push_back( this->value( x,y,z ) );

    tmp2.resize( tmp1.size(),0.0 );

    int i = 0;
    for ( ; i < _iterations; ++i )
    {
        int j = 0;

#pragma omp parallel for
        for ( int z = 0; z < ( int )di_.z; ++z )
            for ( unsigned int y = 0; y < di_.y; ++y )
                for ( unsigned int x = 0; x < di_.x; ++x )
                {
                    if ( i%2==0 ) bilateral_update_simple( tmp1,tmp2,x,y,z,_sigma_r );
                    else        bilateral_update_simple( tmp2,tmp1,x,y,z,_sigma_r );


                    ++j;
                    if ( j % 10000 == 0 )
                        std::cout  << j << " | " << this->n_points() << "                           \r"; std::cout.flush();

                }
        std::cout << "Iteration: " << i << std::endl;
    }

    //store new pixeldata
    for ( unsigned int z = 0; z < di_.z; ++z )
        for ( unsigned int y = 0; y < di_.y; ++y )
            for ( unsigned int x = 0; x < di_.x; ++x )
            {
                if ( i%2==0 )
                {
                    this->value( x,y,z ) = tmp1[di_.x*di_.y*z + di_.x*y + x];
                }
                else
                {
                    this->value( x,y,z ) = tmp2[di_.x*di_.y*z + di_.x*y + x];
                }
            }

    unsigned int elapsed = stop_watch.elapsed();
    std::cout << "Time for filtering: " << elapsed << "ms" << std::endl;
}

//=============================================================================

template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::bilateral_update_simple( std::vector<Scalar> &_src,
        std::vector<Scalar> &_dst,
        int _x, int _y, int _z, Scalar _sigma_r )
{
    int ic = _z*di_.x*di_.y + _y * di_.x + _x;

    Scalar I_center = _src[ic];
    if ( _x >= ( int )di_.x-3 || _x <= 2 ||
            _y >= ( int )di_.y-3  || _y <= 2 ||
            _z >= ( int )di_.z-3  || _z <= 2 ||
            I_center < threshold_ )
    {
        _dst[ic] = 0.0;
        return;
    }

    Scalar acc_w = 0.0;
    Scalar I = 0.0;
    Scalar weight = 0.0;

    Scalar pre_w_int = 1.0/( _sigma_r*sqrt( 2.0 * M_PI ) );
    Scalar quotient =  2.0*_sigma_r*_sigma_r;

    for ( unsigned int iter = 0; iter < cur_gauss_kernel_.size(); ++iter )
    {
        // get precalculated position
        int x = _x + cur_gauss_kernel_[iter].x;
        int y = _y + cur_gauss_kernel_[iter].y;
        int z = _z + cur_gauss_kernel_[iter].z;

        if ( x >= ( int )di_.x || x < 0 ||
                y >= ( int )di_.y || y < 0 ||
                z >= ( int )di_.z || z < 0 ) continue;

        int in = z*di_.x*di_.y + y * di_.x + x;

        // evaluate weight for intensity values
        Scalar diff_int = _src[in] - I_center;
        Scalar w_int = pre_w_int * exp( -( diff_int * diff_int / quotient ) );

        // evaluate product of weights
        weight = w_int * cur_gauss_kernel_[iter].w_dist;

        // increase accumulated weight
        acc_w  += weight;
        I += _src[in]*weight;
    }

    // set new intensity value
    I /= acc_w;
    _dst[ic] = I;
}

//=============================================================================

template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::set_border()
{
    unsigned int x;
    unsigned int y;
    unsigned int z;

    for ( z = 0; z < di_.z; ++z )
    {
        this->value( 0      ,0      ,z ) = 0.0;
        this->value( di_.x-1,0      ,z ) = 0.0;
        this->value( di_.x-1,di_.y-1,z ) = 0.0;
        this->value( 0      ,di_.y-1,z ) = 0.0;
    }

    for ( y = 0; y < di_.y; ++y )
    {
        this->value( 0      ,y      ,0 ) = 0.0;
        this->value( di_.x-1,y      ,di_.z-1 ) = 0.0;
        this->value( di_.x-1,y      ,0 ) = 0.0;
        this->value( 0      ,y      ,di_.z-1 ) = 0.0;
    }

    for ( x = 0; x < di_.x; ++x )
    {
        this->value( x,0      ,0 ) = 0.0;
        this->value( x,di_.y-1,0 ) = 0.0;
        this->value( x,di_.y-1,di_.z-1 ) = 0.0;
        this->value( x,0      ,di_.z-1 ) = 0.0;
    }


}


//=============================================================================


template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::set_bandwidth( Scalar _hs, Scalar _hr )
{
    Vec3 origin( di_.tx, di_.ty, di_.tz );
    Vec3 x_axis( Scalar( di_.x )*di_.sx,0.0,0.0 );
    Vec3 y_axis( 0.0,Scalar( di_.y )*di_.sy,0.0 );
    Vec3 z_axis( 0.0,0.0,Scalar( di_.z )*di_.sz );

    initialize( origin, x_axis/_hs, y_axis/_hs, z_axis/_hs, di_.x, di_.y, di_.z );

    threshold_ = (( max_value_ - min_value_ )*THRESHOLD + min_value_ )/_hr;
    //     std::cout << max_value_ << " , " << min_value_ << " , " << threshold_;
    //     exit(1);

    for ( unsigned int z = 0; z < di_.z; ++z )
        for ( unsigned int y = 0; y < di_.y; ++y )
            for ( unsigned int x = 0; x < di_.x; ++x )
                this->value( x,y,z ) /= _hr;



}

//==============================================================================


template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::init( const Vec3&  _origin,
                            const Vec3&  _x_axis,
                            const Vec3&  _y_axis,
                            const Vec3&  _z_axis,
                            unsigned int            _x_res,
                            unsigned int            _y_res,
                            unsigned int            _z_res )
{
    this->initialize( _origin, _x_axis, _y_axis, _z_axis, _x_res, _y_res, _z_res );
    grad_field_.initialize( _origin, _x_axis, _y_axis, _z_axis, _x_res, _y_res, _z_res );
//     this->resize();
    di_.x = _x_res;
    di_.y = _y_res;
    di_.z = _z_res;
    di_.tx = _origin[0];
    di_.ty = _origin[1];
    di_.tz = _origin[2];
    di_.sx = 1.0;
    di_.sy = 1.0;
    di_.sz = 1.0;
    di_.scale = 1.0;
}

//==============================================================================

template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::normalize()
{
    QTime stop_watch;
    stop_watch.start();

    std::cout << this->values_.size() << std::endl;

    update_min_max();
    normalize( min_value(), max_value() );

    unsigned int elapsed = stop_watch.elapsed();
    std::cout << "Time for normalization: " << elapsed << "ms" << std::endl;
}

//==============================================================================

template< class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::
normalize( double _min, double _max )
{
    unsigned int x_res = di_.x;
    unsigned int y_res = di_.y;
    unsigned int z_res = di_.z;

    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
                this->value( x,y,z ) = ( this->value( x,y,z ) - _min )/( _max - _min );
}

//==============================================================================

template< class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::
clamp( double _min, double _max )
{
    unsigned int x_res = di_.x;
    unsigned int y_res = di_.y;
    unsigned int z_res = di_.z;

    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                if ( this->value( x,y,z ) < _min ) this->value( x,y,z ) = _min;
                if ( this->value( x,y,z ) > _max ) this->value( x,y,z ) = _max;
            }
}

//==============================================================================

template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::getByteDataGradient( std::vector<unsigned char> &_result ) const
{
    unsigned int x_res( this->x_resolution() );
    unsigned int y_res( this->y_resolution() );
    unsigned int z_res( this->z_resolution() );

    Scalar average = 0.0;

    //calculate average of gradient norms
    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
                average+= grad_field_( x,y,z ).norm();

    average /= ( Scalar )x_res*y_res*z_res;

    //calculate standard derivation
    Scalar sigma( 0.0 );
    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
                sigma += ( grad_field_( x,y,z ).norm() - average ) * ( grad_field_( x,y,z ) - average );

    sigma /= ( Scalar )( x_res*y_res*z_res - 1 );
    sigma = sqrt( sigma );

    //clamping data
    Scalar epsilon = 0.0000001;
    Scalar scale = 1.0/( average+1.5*sigma );
    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                if ( grad_field_( x,y,z ) > ( average+1.5*sigma ) && grad_field_( x,y,z ) > epsilon )
                {
                    grad_field_( x,y,z ).normalize();
                }
                else
                {
                    if ( grad_field_( x,y,z ).norm() > epsilon )
                    {
                        grad_field_( x,y,z ).normalize();
                        grad_field_( x,y,z )*=scale;
                    }
                    else grad_field_( x,y,z ) = Vec3( 0.0,0.0,0.0 );
                }
            }

    //normalize all vectors
    //calculate grad_norm
    //for texture: Vec3r(0.0,0.0,0.0) = rgb(128,128,128);
    Scalar grad_norm( 0.0 );
    std::vector<Scalar> grad_norms;
    grad_norms.clear();
    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                grad_norm = grad_field_( x,y,z ).norm();
                if ( grad_norm > epsilon )
                {

                    grad_field_( x,y,z ).normalize();
                    grad_field_( x,y,z ) *= 0.5;
                    grad_field_( x,y,z ) += Vec3( 0.5,0.5,0.5 );
                }
                else
                {
                    grad_field_( x,y,z ) = Vec3( 0.5,0.5,0.5 );
                    grad_norm = 0.0;
                }
                grad_norms.push_back( grad_norm );
            }

    //quantize normal data to unsigned char
    _result.clear();
    _result.resize( 4*x_res*y_res*z_res );

    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                int i = x_res*y_res*z + x_res*y + x;
                Vec3 grad = grad_field_( x,y,z );
                grad_norm     = grad_norms[i];
                _result[4*i  ] = ( unsigned char )( grad[0]   * 255 );
                _result[4*i+1] = ( unsigned char )( grad[1]   * 255 );
                _result[4*i+2] = ( unsigned char )( grad[2]   * 255 );
                _result[4*i+3] = ( unsigned char )( grad_norm * 255 );
            }
}

//==============================================================================

template <class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::getDataGradient( std::vector<float> &_result ) const
{
    unsigned int x_res( this->x_resolution() );
    unsigned int y_res( this->y_resolution() );
    unsigned int z_res( this->z_resolution() );

    _result.clear();
    _result.resize( 4*x_res*y_res*z_res );
    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                // get gradient
                Vec3   grad = grad_field_( x,y,z );
                Scalar norm = grad.norm();

                if ( norm < 1E-6 )
                {
                    grad = Vec3( 0.0,0.0,0.0 );
                    norm = 0.0;
                }
                else
                    grad /= norm;

                int i = x_res*y_res*z + x_res*y + x;
                _result[4*i  ] = ( float )grad[0];
                _result[4*i+1] = ( float )grad[1];
                _result[4*i+2] = ( float )grad[2];
                _result[4*i+3] = ( float )norm;
            }
}

//==============================================================================

template<class Scalar, class Vec3>
Vec3
VolumeImageT<Scalar,Vec3>::calcGradient( const unsigned int &_x,
                                    const unsigned int &_y,
                                    const unsigned int &_z )
{
    if ( _x == 0 || _y == 0 || _z == 0 || _x == this->x_resolution()-1 || _y == this->y_resolution()-1 || _z == this->z_resolution()-1 )
        return Vec3( 0.0,0.0,0.0 );

    Scalar v_000( this->value( _x-1,_y-1,_z-1 ) );
    Scalar v_001( this->value( _x-1,_y-1,_z ) );
    Scalar v_002( this->value( _x-1,_y-1,_z+1 ) );
    Scalar v_010( this->value( _x-1,_y  ,_z-1 ) );
    Scalar v_011( this->value( _x-1,_y  ,_z ) );
    Scalar v_012( this->value( _x-1,_y  ,_z+1 ) );
    Scalar v_020( this->value( _x-1,_y+1,_z-1 ) );
    Scalar v_021( this->value( _x-1,_y+1,_z ) );
    Scalar v_022( this->value( _x-1,_y+1,_z+1 ) );

    Scalar v_100( this->value( _x  ,_y-1,_z-1 ) );
    Scalar v_101( this->value( _x  ,_y-1,_z ) );
    Scalar v_102( this->value( _x  ,_y-1,_z+1 ) );
    Scalar v_110( this->value( _x  ,_y  ,_z-1 ) );
    Scalar v_112( this->value( _x  ,_y  ,_z+1 ) );
    Scalar v_120( this->value( _x  ,_y+1,_z-1 ) );
    Scalar v_121( this->value( _x  ,_y+1,_z ) );
    Scalar v_122( this->value( _x  ,_y+1,_z+1 ) );

    Scalar v_200( this->value( _x+1,_y-1,_z-1 ) );
    Scalar v_201( this->value( _x+1,_y-1,_z ) );
    Scalar v_202( this->value( _x+1,_y-1,_z+1 ) );
    Scalar v_210( this->value( _x+1,_y  ,_z-1 ) );
    Scalar v_211( this->value( _x+1,_y  ,_z ) );
    Scalar v_212( this->value( _x+1,_y  ,_z+1 ) );
    Scalar v_220( this->value( _x+1,_y+1,_z-1 ) );
    Scalar v_221( this->value( _x+1,_y+1,_z ) );
    Scalar v_222( this->value( _x+1,_y+1,_z+1 ) );

    Scalar w1 = 1.0;
    Scalar w2 = 2.0;
    Scalar w3 = 4.0;


    Scalar dfdx = ( w1*v_200 + w2*v_201 + w1*v_202
                    + w2*v_210 + w3*v_211 + w2*v_212
                    + w1*v_220 + w2*v_221 + w1*v_222 )
                  -( w1*v_000 + w2*v_001 + w1*v_002
                     + w2*v_010 + w3*v_011 + w2*v_012
                     + w1*v_020 + w2*v_021 + w1*v_022 );

    Scalar dfdy = ( w1*v_020 + w2*v_120 + w1*v_220
                    + w2*v_021 + w3*v_121 + w2*v_221
                    + w1*v_022 + w2*v_122 + w1*v_222 )
                  -( w1*v_000 + w2*v_100 + w1*v_200
                     + w2*v_001 + w3*v_101 + w2*v_201
                     + w1*v_002 + w2*v_102 + w1*v_202 );

    Scalar dfdz = ( w1*v_002 + w2*v_102 + w1*v_202
                    + w2*v_012 + w3*v_112 + w2*v_212
                    + w1*v_022 + w2*v_122 + w1*v_222 )
                  -( w1*v_000 + w2*v_100 + w1*v_200
                     + w2*v_010 + w3*v_110 + w2*v_210
                     + w1*v_020 + w2*v_120 + w1*v_220 );

    dfdx /= di_.sx * ( 8.0 * w1 + 8.0* w2 + 2.0*w3 );
    dfdy /= di_.sy * ( 8.0 * w1 + 8.0* w2 + 2.0*w3 );
    dfdz /= di_.sz * ( 8.0 * w1 + 8.0* w2 + 2.0*w3 );


    return Vec3( dfdx, dfdy, dfdz );
}

//==============================================================================


template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::assign_gradients()
{
    unsigned int x_res( this->x_resolution() );
    unsigned int y_res( this->y_resolution() );
    unsigned int z_res( this->z_resolution() );

#pragma omp parallel for
    for ( int z = 0; z < ( int )z_res; ++z )
    {
        std::cout << Scalar( z ) / Scalar( z_res ) * 100.0 << "% ready...\r";
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
//                 int idx = z * x_res * y_res + y * x_res + x;
                grad_field_( x,y,z ) = -calcGradient( x,y,z );;
            }
    }
    std::cout << "100.0 % ready...\n";
}

//==============================================================================


template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::normalize_gradients()
{
    unsigned int x_res( this->x_resolution() );
    unsigned int y_res( this->y_resolution() );
    unsigned int z_res( this->z_resolution() );

    Scalar min = std::numeric_limits<Scalar>::max();
    Scalar max = 0.0;

    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                Scalar grad = grad_field_( x,y,z ).norm();
                if ( grad > max ) max = grad;
                else if ( grad < min ) min = grad;
            }

    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                Vec3 grad = grad_field_( x,y,z );
                Scalar scale = ( grad.norm() - min )/( max - min );
                if ( grad.norm() > 0.000001 )
                {
                    grad.normalize();
                    grad *= scale;
                    grad_field_( x,y,z ) = grad;
                }
                else grad_field_( x,y,z ) = Vec3( 0.0,0.0,0.0 );
            }
}

//==============================================================================


template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::histogram( const int _n_bins, std::vector<unsigned int> &_bins, Scalar &_size, bool _to_image )
{
    unsigned int x_res( this->x_resolution() );
    unsigned int y_res( this->y_resolution() );
    unsigned int z_res( this->z_resolution() );

    std::vector<unsigned int>::iterator b_iter;

    _bins.clear();
    _bins.resize( _n_bins,0 );

    Scalar imax = -std::numeric_limits<Scalar>::max();
    Scalar imin = std::numeric_limits<Scalar>::max();
    Scalar val;

    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                val = Scalar( this->value( x,y,z ) );
                if ( val > imax ) imax = val;
                else if ( val < imin ) imin = val;
            }

    _size = ( imax-imin )/Scalar( _n_bins );

    unsigned int bin_index = 0;
    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                val = this->value( x,y,z );
                bin_index = ( unsigned int )(( val-imin )/( imax-imin )*Scalar( _n_bins-1 ) );
                ++_bins.at( bin_index );
            }

    ///Image of Histogram
    if ( _to_image )
    {
        QImage histogram( 5*_n_bins, 5*_n_bins, QImage::Format_RGB32 );
        histogram.fill( 0 );

        unsigned int bmax = 0;
        unsigned int bmin = std::numeric_limits<unsigned int>::max();
        b_iter = _bins.begin();
        for ( ; b_iter != _bins.end(); ++b_iter )
        {
            if (( *b_iter ) > bmax ) bmax = ( *b_iter );
            if (( *b_iter ) < bmin ) bmin = ( *b_iter );
        }

        std::cout << "Histogram" << std::endl;
        std::cout << "---------" << std::endl;
        std::cout << "Min Intensity: " << std::setw( 10 ) << imin << ",   ";
        std::cout << "Max Intensity: " << std::setw( 10 ) << imax << std::endl;
        std::cout << "Min Entries  : " << std::setw( 10 ) << bmin << ",   ";
        std::cout << "Max Entries  : " << std::setw( 10 ) << bmax << std::endl << std::endl;

        int x = 0;
        int y = 0;
        Scalar scale = Scalar( _n_bins-1 )*5.0;
        for ( unsigned int bi = 0; bi < _bins.size(); ++bi )
        {
            if ( _bins.at( bi ) > x_res*y_res*z_res*0.0001 )
            {
                x = 5*bi;
                y = int( scale*Scalar( _bins.at( bi ) - bmin )/Scalar( bmax-bmin ) );
                for ( int dy = -2; dy <= 2; ++dy )
                {
                    if (( y+dy ) >= 0 && ( y+dy ) < 5*_n_bins )
                    {
                        histogram.setPixel( x,( 5*_n_bins-1 )-( y+dy ),qRgb( 255,255,255 ) );
                    }
                }
                for ( int dx = -2; dx <= 2; ++dx )
                {
                    if (( x+dx ) >= 0 && ( x+dx ) < 5*_n_bins )
                    {
                        histogram.setPixel( x+dx,( 5*_n_bins-1 )-y,qRgb( 255,255,255 ) );
                    }
                }
            }
        }
        histogram.save( "Histogram.png","PNG" );
    }
}

//==============================================================================


template<class Scalar, class Vec3>
void
VolumeImageT<Scalar,Vec3>::normalize_with_histogram( Scalar _max, Scalar _percent, int _n_bins )
{
    unsigned int x_res( this->x_resolution() );
    unsigned int y_res( this->y_resolution() );
    unsigned int z_res( this->z_resolution() );

    Scalar size = 0.0;
    std::vector<unsigned int> bins;

    Scalar imax = -std::numeric_limits<Scalar>::max();
    Scalar imin = std::numeric_limits<Scalar>::max();
    Scalar val;


    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                val = Scalar( this->value( x,y,z ) );
                if ( val > imax ) imax = val;
                else if ( val < imin ) imin = val;
            }

    histogram( _n_bins,bins,size );

    unsigned int sum = 0;
    for ( unsigned int i = 0; i < bins.size(); ++i )
        sum += bins[i];

    std::cout << sum << std::endl;

    sum = 0;

    bool sum_low = true;
    unsigned int b_iter = 0;
    while ( sum_low && b_iter < bins.size() )
    {
        sum += bins.at( b_iter );
        if ( sum < ( unsigned int )Scalar( x_res*y_res*z_res )*_percent )
        {
            sum_low = true;
            imin = Scalar( Scalar( b_iter )*size );
        }
        else sum_low = false;
        ++b_iter;
    }
    std::cout << b_iter << std::endl << std::endl;

    sum = 0;
    sum_low = true;
    b_iter = bins.size()-1 ;
    while ( sum_low && b_iter >= 0 )
    {
        sum += bins.at( b_iter );
        if ( sum < ( unsigned int )Scalar( x_res*y_res*z_res )*_percent )
        {
            sum_low = true;
            imax = Scalar( Scalar( b_iter )*size );
        }
        else sum_low = false;
        --b_iter;
    }

    std::cout << b_iter << std::endl << std::endl;
    std::cout << imin << ", " << imax << std::endl;

    for ( unsigned int z = 0; z < z_res; ++z )
        for ( unsigned int y = 0; y < y_res; ++y )
            for ( unsigned int x = 0; x < x_res; ++x )
            {
                val = this->value( x,y,z );

                if ( val > imax ) val = imax;
                else if ( val < imin ) val = imin;
                else
                {
                    //set data
                    val = Scalar( _max*Scalar( val - imin )/Scalar( imax - imin ) );
                    this->value( x,y,z ) = val;
                }
            }
}

//==============================================================================
} // end namespace
//==============================================================================

//==============================================================================

// Local Variables:
// mode: C++
// End:

//==============================================================================


