/*
* This code is in the public domain.
* 
* Author: Christopher Gary
* 
* All credit for algorithm development belongs to the authors of LAPACK.
* */

#include <cmath>
#include <cctype>
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <utility>
#include <vector>
#include <limits>

using std::min;
using std::max;
using std::isnan;

// Translations with index-space refactoring from Reference-LAPACK.
// 
// Refer to the original FORTRAN routines for documentation.
//
// Several routines were changed to facilitate a standalone implementation.
//
// This is intended for reference purposes only; not meant for production code.
//
// The only guarantee here is that everything is 0-based, and should be
// easier to follow if writing -something new- for a different value type or
// something generic around the same algorithms.
struct Toy_Cxx_LAPACK_3_7_0
{
  struct XerblaException
  {
    const char *name;
    int info;
  };

  static constexpr int col_mjr( int i, int j, int ld ) noexcept
  { return i + j*ld; }

  // This may or may not have improved readability...
  static constexpr inline double sqr( double x ) noexcept
  { return x*x; }

  static /*constexpr*/ inline bool dlaeisnan( double u, double v ) noexcept
  { return std::isnan(u) || std::isnan(v); /*u != v;*/ }

  static /*constexpr*/ inline double dlapy2( double x, double y ) noexcept
  { return hypot( x, y ); }

  enum TransposeOpt : char
  {
    kTrans = 'T',
    kNoTrans = 'N',
    kCeeTrans_idk = 'C'
  };

  enum SideOpt : char
  {
    kLeft = 'L',
    kRight = 'R'
  };

  enum HalfOpt : char
  {
    kUpper = 'U',
    kLower = 'L',
    kUpAndLo = 'F'
  };

  enum TypeOpt : char
  {
    kFullMat = 'G',
    kLowerTri = 'L',
    kUpperTri = 'U',
    kUpperHess = 'H',
    kLowerSyBand = 'B',
    kUpperSyBand = 'Q',
    kBandMat = 'Z'
  };

  enum NormOpt : char
  {
    kMaxNorm = 'M',
    kOneNorm = '1',
    kInfNorm = 'I',
    kFrobNorm = 'F'
  };

  enum DiagIsUnitOpt : char
  {
    kDiagUnit = 'U',
    kDiagNotUnit = 'N'
  };

  enum PivotOpt : char
  {
    kPivotBtm = 'B',
    kPivotTop = 'T',
    kPivotVar = 'V'
  };

  enum DirectOpt : char
  {
    kForwd = 'F',
    kBackwd = 'B'
  };

  enum StoreOpt : char
  {
    kByCol = 'C',
    kByRow = 'R'
  };

  enum CompEigvOpt : char
  {
    kNoEigv = 'N',  // Eigenvectors only
    kOrigEigv = 'V',// Eigensystem of the orignal matrix
    kTridEigv = 'I' // Eigensystem of the tridiagonal matrix
  };

  enum SortDirectOpt : char
  {
    kSortInc = 'I',
    kSortDec = 'D'
  };


  //------------------------------------------------------------------------
  // dcopy
  //------------------------------------------------------------------------

  void dcopy( int n, double *dx, int incx, double *dy, int incy )
  {
    int i, ix, iy;

    if( n <= 0 )
    { return; }

    if( ( 1 == incx ) && ( 1 == incy ) )
    {
      memcpy( dy, dx, n*sizeof(double) );
    }
    else
    {
      ix = 0;
      iy = 0;
      if( incx < 0 ){ ix = ((-n+1)*incx+1)-1; }
      if( incy < 0 ){ iy = ((-n+1)*incy+1)-1; }
      for( i = 0; i < n; ++i )
      {
        dy[iy] = dx[ix];
        ix += incx;
        iy += incy;
      }
    }
  }

  //------------------------------------------------------------------------
  // dcombssq
  //------------------------------------------------------------------------

  void dcombssq( double v1[2], double v2[2] )
  {
    if( v1[0] >= v2[0] )
    {
      if( 0.0 != v1[0] )
      { v1[1] = v2[1] + sqr(v2[0]/v1[0])*v2[1]; }else
      { v1[1] = v1[1] + v2[1]; }
    }
    else
    {
      v1[1] = v2[1] + sqr(v1[0]/v2[0])*v1[1];
      v1[0] = v2[0];
    }
  }

  //------------------------------------------------------------------------
  // dlassq
  //------------------------------------------------------------------------

  void dlassq( int n, double *x, int incx, double &scale, double &sumsq )
  {
    int ix;
    double absxi;

    if( n > 0 )
    {
      for( ix = 0; ix < (1 + (n-1)*incx); ix += incx )
      {
        absxi = fabs( x[ix] );
        if( ( absxi > 0.0 ) || isnan(absxi) )
        {
          if( scale < absxi )
          {
            sumsq = 1.0 + sumsq*sqr( scale / absxi );
            scale = absxi;
          }
          else
          {
            sumsq += sqr( absxi / scale );
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // daxpy
  //------------------------------------------------------------------------

  void daxpy( int n, double da, double *dx, int incx, double *dy, int incy )
  {
    int i, ix, iy;

    if( n <= 0 ){ return; }
    if( 0.0 == da ){ return; }
    if( ( 1 == incx ) && ( 1 == incy ) )
    {
      // clang should auto-vectorize things like this with -O3...
      for( i = 0; i < n; ++i )
      { dy[i] += da*dx[i]; }
    }
    else
    {
      ix = 0;
      iy = 0;
      if( incx < 0 ){ ix = (((-n+1)*incx)+1)-1; }
      if( incy < 0 ){ iy = (((-n+1)*incy)+1)-1; }
      for( i = 0; i < n; ++i )
      {
        dy[iy] = da*dx[ix];
        ix += incx;
        iy += incy;
      }
    }
  }

  //------------------------------------------------------------------------
  // ddot
  //------------------------------------------------------------------------

  double ddot( int n, double *dx, int incx, double *dy, int incy )
  {
    double dtmp;
    int i, ix, iy;

    if( n <= 0 ){ return 0.0; }

    dtmp = 0.0;

    if( ( 1 == incx ) && ( 1 == incy ) )
    {
      for( i = 0; i < n; ++i )
      { dtmp += dx[i]*dy[i]; }
    }
    else
    {
      ix = 0;
      iy = 0;
      if( incx < 0 ){ ix = (((-n+1)*incx)+1)-1; }
      if( incy < 0 ){ iy = (((-n+1)*incy)+1)-1; }
      for( i = 0; i < n; ++i )
      {
        dtmp += dx[ix]*dy[iy];
        ix += incx;
        iy += incy;
      }
    }

    return dtmp;
  }

  //------------------------------------------------------------------------
  // dscal
  //------------------------------------------------------------------------

  void dscal( int n, double da, double *dx, int incx )
  {
    int i, nincx;

    if( (n <= 0) || (incx <= 0) ){ return; }

    nincx = n*incx;
    for( i = 0; i < nincx; i += incx )
    { dx[i] = da*dx[i]; }
  }

  //------------------------------------------------------------------------
  // dlae2
  //------------------------------------------------------------------------

  void dlae2( double a, double b, double c, double &rt1, double &rt2 )
  {
    double ab, acmn, acmx, adf, df, rt, sm, tb;

    // Compute the eigenvalues of [ a b ]
    //                            [ b c ]

    sm = a + c;
    df = a - c;
    adf = fabs( df );
    tb = b + b;
    ab = fabs( tb );

    if( fabs(a) > fabs(c) )
    {
      acmx = a;
      acmn = c;
    }
    else
    {
      acmx = c;
      acmn = a;
    }

    if( adf > ab )
    {
      rt = adf*sqrt( 1.0 + sqr(ab/adf) );
    }
    else if( adf < ab )
    {
      rt = ab*sqrt( 1.0 + sqr(adf/ab) );
    }
    else
    {
      // Includes case AB=ADF=0

      rt = ab*sqrt(2.0);
    }

    if( sm < 0.0 )
    {
      rt1 = 0.5*( sm - rt );
      rt2 = (acmx/rt1)*acmn - (b/rt1)*b;
    }
    else if( sm > 0.0 )
    {
      rt1 = 0.5*( sm + rt );
      rt2 = (acmx/rt1)*acmn - (b/rt1)*b;
    }
    else
    {
      // Includes case RT1 = RT2 = 0
      rt1 = 0.5*rt;
      rt2 = -0.5*rt;
    }
  }

  //------------------------------------------------------------------------
  // dlaev2
  //------------------------------------------------------------------------

  void dlaev2( double a, double b, double c, double &rt1, double &rt2, double &cs1, double &sn1 )
  {
    int sgn1, sgn2;
    double ab, acmn, acmx, acs, adf, cs, ct, df, rt, sm, tb, tn;

    // Compute the eigenvalues

    sm = a + c;
    df = a - c;
    adf = fabs( df );
    tb = b + b;
    ab = fabs( tb );

    if( fabs(a) > fabs(c) )
    {
      acmx = a;
      acmn = c;
    }
    else
    {
      acmx = c;
      acmn = a;
    }

    if( adf > ab )
    {
      rt = adf*sqrt( 1.0 + sqr(ab/adf) );
    }
    else if( adf < ab )
    {
      rt = ab*sqrt( 1.0 + sqr(adf/ab) );
    }
    else
    {
      // Includes case AB=ADF=0

      rt = ab*sqrt(2.0);
    }

    if( sm < 0.0 )
    {
      rt1 = 0.5*( sm - rt );
      sgn1 = -1;
      rt2 = (acmx/rt1)*acmn - (b/rt1)*b;
    }
    else if( sm > 0.0 )
    {
      rt1 = 0.5*( sm + rt );
      sgn1 = 1;
      rt2 = (acmx/rt1)*acmn - (b/rt1)*b;
    }
    else
    {
      // Includes case RT1 = RT2 = 0
      rt1 = 0.5*rt;
      rt2 = -0.5*rt;
      sgn1 = 1;
    }

    // Compute the eigenvector

    if( df >= 0.0 )
    {
      cs = df + rt;
      sgn2 = 1;
    }
    else
    {
      cs = df - rt;
      sgn2 = -1;
    }

    acs = fabs(cs);

    if( acs > ab )
    {
      ct = -tb/cs;
      sn1 = 1.0/sqrt( 1.0 + sqr(ct) );
      cs1 = ct*sn1;
    }
    else
    {
      if( 0.0 == ab )
      {
        cs1 = 1.0;
        sn1 = 0.0;
      }
      else
      {
        tn = -cs/tb;
        cs1 = 1.0/sqrt( 1.0 + sqr(tn) );
        sn1 = tn*cs1;
      }
    }
    if( sgn1 == sgn2 )
    {
      tn = cs1;
      cs1 = -sn1;
      sn1 = tn;
    }
  }

  //------------------------------------------------------------------------
  // dlartg
  //------------------------------------------------------------------------

  void dlartg( double f, double g, double &cs, double &sn, double &r )
  {
    int count;
    double f1, g1, safmin, safmax, eps, scale;

    // NOTE: Probably slightly wrong...
    safmin = DBL_MIN;
    safmax = DBL_MAX;
    eps = DBL_EPSILON;

    if( 0.0 == g )
    {
      cs = 1.0;
      sn = 0.0;
      r = f;
    }
    else if( 0.0 == f )
    {
      cs = 0.0;
      sn = 1.0;
      r = g;
    }
    else
    {
      f1 = f;
      g1 = g;
      scale = max( fabs(f1), fabs(g1) );
      if( scale >= safmax )
      {
        count = 0;
        do
        {
          ++count;
          f1 *= safmin;
          g1 *= safmin;
          scale = max( fabs(f1), fabs(g1) );
        }
        while( safmax >= scale );

        r = hypot( f1, g1 );
        cs = f1/r;
        sn = g1/r;
        while( count-- )
        { r *= safmax; }
      }
      else if( scale <= safmin )
      {
        count = 0;
        do
        {
          ++count;
          f1 *= safmax;
          g1 *= safmax;
          scale = max( fabs(f1), fabs(g1) );
        }
        while( scale <= safmin );

        r = hypot( f1, g1 );
        cs = f1/r;
        sn = g1/r;
        while( count-- )
        { r *= safmin; }
      }
      else
      {
        r = hypot( f1, g1 );
        cs = f1/r;
        sn = g1/r;
      }

      if( (fabs(f) > fabs(g)) && (cs < 0.0) )
      {
        cs = -cs;
        sn = -sn;
        r = -r;
      }
    }
  }

  //------------------------------------------------------------------------
  // dlaset
  //------------------------------------------------------------------------

  void dlaset( HalfOpt uplo, int m, int n, double alpha, double beta, double *pA, int lda )
  {
    int i, j;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    if( kUpper == uplo )
    {
      // Set the strictly upper triangular or trapezoidal part of the
      // array to ALPHA.
      for( j = 1; j < n; ++j )
      { for( i = 0; i < min(j-1,m-1); ++i )
      { A(i,j) = alpha; } }
    }
    else if( kLower == uplo )
    {
      // Set the strictly lower triangular or trapezoidal part of the
      // array to ALPHA.
      for( j = 0; j < min(m,n); ++j )
      { for( i = j+1; i < m; ++i )
      { A(i,j) = alpha; } }
    }
    else
    {
      // Set the leading m-by-n submatrix to ALPHA.
      for( j = 0; j < n; ++j )
      { for( i = 0; i < m; ++i )
      { A(i,j) = alpha; } }
    }

    // Set the first min(M,N) diagonal elements to BETA.
    for( i = 0; i < min(m,n); ++i )
    { A(i,i) = beta; }
  }

  //------------------------------------------------------------------------
  // dlasr
  //------------------------------------------------------------------------

  void dlasr( SideOpt side, PivotOpt pivot, DirectOpt direct, int m, int n, double *c, double *s, double *pA, int lda )
  {
    int i, j;
    double ctmp, stmp, tmp;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Test the input parameters

    if( m < 0 ){ throw XerblaException{ "dlasr", 4 }; }
    if( n < 0 ){ throw XerblaException{ "dlasr", 5 }; }
    if( lda < max(1,m) ){ throw XerblaException{ "dlasr", 9 }; }

    // Quick return if possible

    if( 0 == m || 0 == n ){ return; }

    if( kLeft == side )
    {
      // Form  P * A

      switch( pivot )
      {
      default: break;
      case kPivotVar:
        {
          if( kForwd == direct )
          {
            for( j = 0; j < (m-1); ++j )
            {
              ctmp = c[j];
              stmp = s[j];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < n; ++i )
                {
                  tmp = A(j+1,i);
                  A(j+1,i) = ctmp*tmp - stmp*A(j,i);
                  A(j,i) = stmp*tmp + ctmp*A(j,i);
                }
              }
            }
          }
          else if( kBackwd == direct )
          {
            for( j = (m-1)-1; j >= 0; --j )
            {
              ctmp = c[j];
              stmp = s[j];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < n; ++i )
                {
                  tmp = A(j+1,i);
                  A(j+1,i) = ctmp*tmp - stmp*A(j,i);
                  A(j,i) = stmp*tmp + ctmp*A(j,i);
                }
              }
            }
          }
        }
        break;

      case kPivotTop:
        {
          if( kForwd == direct )
          {
            for( j = 1; j < m; ++j )
            {
              ctmp = c[j-1];
              stmp = s[j-1];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < n; ++i )
                {
                  tmp = A(j,i);
                  A(j,i) = ctmp*tmp - stmp*A(0,i);
                  A(0,i) = stmp*tmp + ctmp*A(0,i);
                }
              }
            }
          }
          else if( kBackwd == direct )
          {
            for( j = m-1; j >= 1; --j )
            {
              ctmp = c[j-1];
              stmp = s[j-1];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < n; ++i )
                {
                  tmp = A(j,i);
                  A(j,i) = ctmp*tmp - stmp*A(0,i);
                  A(0,i) = stmp*tmp + ctmp*A(0,i);
                }
              }
            }
          }
        }
        break;

      case kPivotBtm:
        {
            if( kForwd == direct )
            {
              for( j = 0; j < m-1; ++j )
              {
                ctmp = c[j];
                stmp = s[j];
                if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
                {
                  for( i = 0; i < n; ++i )
                  {
                    tmp = A(j,i);
                    A(j,i) = stmp*A(m-1,i) + ctmp*tmp;
                    A(m-1,i) = ctmp*A(m-1,i) - stmp*tmp;
                  }
                }
              }
            }
            else if( kBackwd == direct )
            {
              for( j = (m-1)-1; j >= 0; --j )
              {
                ctmp = c[j];
                stmp = s[j];
                if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
                {
                  for( i = 0; i < n; ++i )
                  {
                    tmp = A(j,i);
                    A(j,i) = stmp*A(m-1,i) + ctmp*tmp;
                    A(m-1,i) = ctmp*A(m-1,i) - stmp*tmp;
                  }
                }
              }
            }
        }
        break;
      }
    }
    else if( kRight == side )
    {
      // Form A * (~P)

      switch( pivot )
      {
      default: break;
      case kPivotVar:
        {
          if( kForwd == direct )
          {
            for( j = 0; j < n-1; ++j )
            {
              ctmp = c[j];
              stmp = s[j];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < m; ++i )
                {
                  tmp = A(i,j+1);
                  A(i,j+1) = ctmp*tmp - stmp*A(i,j);
                  A(i,j) = stmp*tmp + ctmp*A(i,j);
                }
              }
            }
          }
          else if( kBackwd == direct )
          {
            for( j = (n-1)-1; j >= 0; --j )
            {
              ctmp = c[j];
              stmp = s[j];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < m; ++i )
                {
                  tmp = A(i,j+1);
                  A(i,j+1) = ctmp*tmp - stmp*A(i,j);
                  A(i,j) = stmp*tmp + ctmp*A(i,j);
                }
              }
            }
          }
        }
        break;

      case kPivotTop:
        {
          if( kForwd == direct )
          {
            for( j = 1; j < n; ++j )
            {
              ctmp = c[j-1];
              stmp = s[j-1];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < m; ++i )
                {
                  tmp = A(i,j);
                  A(i,j) = ctmp*tmp - stmp*A(i,0);
                  A(i,0) = stmp*tmp + ctmp*A(i,0);
                }
              }
            }
          }
          else if( kBackwd == direct )
          {
            for( j = n-1; j >= 1; --j )
            {
              ctmp = c[j-1];
              stmp = s[j-1];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < m; ++i )
                {
                  tmp = A(i,j);
                  A(i,j) = ctmp*tmp - stmp*A(i,0);
                  A(i,0) = stmp*tmp + ctmp*A(i,0);
                }
              }
            }
          }
        }
        break;

      case kPivotBtm:
        {
          if( kForwd == direct )
          {
            for( j = 0; j < n-1; ++j )
            {
              ctmp = c[j];
              stmp = s[j];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < m; ++i )
                {
                  tmp = A(i,j);
                  A(i,j) = stmp*A(i,n-1) + ctmp*tmp;
                  A(i,n-1) = ctmp*A(i,n) - stmp*tmp;
                }
              }
            }
          }
          else if( kBackwd == direct )
          {
            for( j = (n-1)-1; j >= 0; --j )
            {
              ctmp = c[j];
              stmp = s[j];
              if( ( 1.0 != ctmp ) || ( 0.0 != stmp ) )
              {
                for( i = 0; i < m; ++i )
                {
                  tmp = A(i,j);
                  A(i,j) = stmp*A(i,n-1) + ctmp*tmp;
                  A(i,n-1) = ctmp*A(i,n-1) - stmp*tmp;
                }
              }
            }
          }
        }
        break;
      }
    }
  }

  //------------------------------------------------------------------------
  // dlasrt
  //------------------------------------------------------------------------

  void dlasrt( SortDirectOpt id, int n, double *d )
  {
    if( n < 0 ){ throw XerblaException{ "dlasrt", 2 }; }

    if( kSortInc == id )
    {
      qsort( d, (size_t)n, sizeof(double),
      []( const void *pa, const void *pb ) noexcept -> int
      {
        const auto &a = *(double*)pa;
        const auto &b = *(double*)pb;
        if( a > b ){ return  1; }
        if( a < b ){ return -1; }
        return 0;
      });
    }
    else
    {
      qsort( d, (size_t)n, sizeof(double),
      []( const void *pa, const void *pb ) noexcept -> int
      {
        const auto &a = *(double*)pa;
        const auto &b = *(double*)pb;
        if( a > b ){ return -1; }
        if( a < b ){ return  1; }
        return 0;
      });
    }
  }

  //------------------------------------------------------------------------
  // iladlc
  //------------------------------------------------------------------------

  int iladlc( int m, int n, double *pA, int lda )
  {
    int i, j;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Quick test for the common case where one corner is non-zero.
    if( 0 == n )
    { return n-1; }
    if( (0.0 != A(0,n-1)) || (0.0 != A(m-1,n-1)) )
    { return n-1; }

    // Now scan each column from the end, returning with the first non-zero.
    for( j = n-1; j >= 0; --j )
    { for( i = 0; i < m; ++i )
    { if( 0.0 != A(i,j) ){ break; } } }

    // NOTE: The original function did not initialize its return value,
    // and would produce an uknown value if the matrix was all zeros.
    //
    // It is assumed all client code is using the function under
    // the correct assumptions.
    return j;
  }

  //------------------------------------------------------------------------
  // iladlr
  //------------------------------------------------------------------------

  int iladlr( int m, int n, double *pA, int lda )
  {
    int i, j, row;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Quick test for the common case where one corner is non-zero.
    if( 0 == m )
    { return m-1; }
    if( (0.0 != A(m-1,0)) || (0.0 != A(m-1,n-1)) )
    { return m-1; }

    row = -1;

    // Scan up each column tracking the last zero row seen.
    for( j = 0; j < n; ++j )
    {
      i = m-1;
      while( (0.0 == A(max(i,0),j)) && (i >= 0) ){ --i; }
      row = max( row, i );
    }

    return row;
  }

  //------------------------------------------------------------------------
  // dlascl
  //------------------------------------------------------------------------

  void dlascl( TypeOpt type, int kl, int ku, double cfrom, double cto, int m, int n, double *pA, int lda )
  {
    bool done;
    int i, j, k1, k2, k3, k4;
    double bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Test the input arguments

    if( (0.0 == cfrom) || isnan(cfrom) ){ throw XerblaException{ "dlascl", 4 }; }
    if( isnan(cto) ){ throw XerblaException{ "dlascl", 5 }; }
    if( m < 0 ){ throw XerblaException{ "dlascl", 6 }; }
    if( (n < 0) ||
      ((kLowerSyBand == type) && (n != m)) ||
      ((kUpperSyBand == type) && (n != m)) ){ throw XerblaException{ "dlascl", 7 }; }
    if( ( type < kUpperHess ) && (lda < max(1,m)) ){ throw XerblaException{ "dlascl", 9 }; }
    if( type >= kLowerSyBand )
    {
      if( (kl < 0) || (kl > max(m-1,0)) ){ throw XerblaException{ "dlascl", 2 }; }

      if( (ku < 0) || (ku > max(n-1,0)) ||
        ((( kLowerSyBand == type ) ||
          ( kUpperSyBand == type )) &&
        (kl != ku)) ){ throw XerblaException{ "dlascl", 3 }; }

      if( (( kLowerSyBand == type ) && (lda < (kl+1))) ||
          (( kUpperSyBand == type ) && (lda < (ku+1))) ||
          (( kBandMat == type ) &&
          (lda < 2*(kl+ku+1))) ){ throw XerblaException{ "dlascl", 9 }; }
    }

    // Quick return if possible

    if( (0 == n) || (0 == m )){ return; }

    smlnum = DBL_MIN;
    bignum = 1.0/smlnum;

    cfromc = cfrom;
    ctoc = cto;

    done = false;

    do
    {
      cfrom1 = cfromc*smlnum;
      if( cfrom1 == cfromc )
      {
        // CFROMC is an inf.  Multiply by a correctly signed zero for
        // finite CTOC, or a NaN if CTOC is infinite.
        mul = ctoc / cfromc;
        done = true;
        cto1 = ctoc;
      }
      else
      {
        cto1 = ctoc / bignum;
        if( cto1 == ctoc )
        {
          // CTOC is either 0 or an inf.  In both cases, CTOC itself
          // serves as the correct multiplication factor.
          mul = ctoc;
          done = true;
          cfromc = 1.0;
        }
        else if( (fabs(cfrom1) > fabs(ctoc)) && (0.0 != ctoc) )
        {
          mul = smlnum;
          done = false;
          cfromc = cfrom1;
        }
        else if( fabs(cto1) > fabs(cfromc) )
        {
          mul = bignum;
          done = false;
          ctoc = cto1;
        }
        else
        {
          mul = ctoc / cfromc;
          done = true;
        }
      }

      switch( type )
      {
      // Technically should not get here.
      default: throw XerblaException{ "dlascl", 1 };

      case kFullMat:
        {
          for( j = 0; j < n; ++j )
          { for( i = 0; i < m; ++i )
          { A(i,j) *= mul; } }
        }
        break;

      case kLowerTri:
        {
          for( j = 0; j < n; ++j )
          { for( i = j; i < m; ++i )
          { A(i,j) *= mul; } }
        }
        break;

      case kUpperTri:
        {
          for( j = 0; j < n; ++j )
          { for( i = 0; i <= min(j,m-1); ++i )
          { A(i,j) *= mul; } }
        }
        break;

      case kUpperHess:
        {
          for( j = 0; j < n; ++j )
          { for( i = 0; i <= min(j+1,m-1); ++i )
          { A(i,j) *= mul; } }
        }
        break;

      case kLowerSyBand:
        {
          k3 = kl + 1;
          k4 = n + 1;
          for( j = 0; j < n; ++j )
          { for( i = 0; i <= (min(k3,k4-(j+1))-1); ++i )
          { A(i,j) *= mul; } }
        }
        break;

      case kUpperSyBand:
        {
          k1 = ku + 2;
          k3 = ku + 1;
          for( j = 0; j < n; ++j )
          { for( i = (max(k1-(j+1),1)-1); i < k3; ++i )
          { A(i,j) *= mul; } }
        }
        break;

      case kBandMat:
        {
          k1 = kl + ku + 2;
          k2 = kl + 1;
          k3 = 2*kl + ku + 1;
          k4 = kl + ku + 1 + m;
          for( j = 0; j < n; ++j )
          { for( i = (max(k1-(j+1),k2)-1); i < min(k3,k4-(j+1)); ++i )
          { A(i,j) *= mul; } }
        }
        break;
      }
    }
    while( ! done );
  }

  //------------------------------------------------------------------------
  // dsyr2
  //------------------------------------------------------------------------

  void dsyr2( HalfOpt uplo, int n, double alpha, double *x, int incx, double *y, int incy, double *pA, int lda )
  {
    double tmp1, tmp2;
    int i, ix, iy, j, jx, jy, kx, ky;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Test the input parameters.

    if( n < 0 ){ throw XerblaException{ "dsyr2", 2 }; }
    if( 0 == incx ){ throw XerblaException{ "dsyr2", 5 }; }
    if( 0 == incy ){ throw XerblaException{ "dsyr2", 7 }; }
    if( lda < max(1,n) ){ throw XerblaException{ "dsyr2", 9 }; }

    // Quick return if possible.

    if( (0 == n) ||
      (0.0 == alpha) ){ return; }

    // Set up the start points in X and Y if the increments are not both
    // unity.

    if( (1 != incx) || (1 != incy) )
    {
      if( incx > 0 )
      { kx = 0; }else
      { kx = (1 - (n-1)*incx) - 1; }
      if( incy > 0 )
      { ky = 0; }else
      { ky = (1 - (n-1)*incy) - 1; }
      jx = kx;
      jy = ky;
    }

    // Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through the triangular part
    // of A.

    if( kUpper == uplo )
    {
      // Form  A  when A is stored in the upper triangle.

      if( (1 == incx) && (1 == incy) )
      {
        for( j = 0; j < n; ++j )
        {
          if( (0.0 != x[j]) || (0.0 != y[j]) )
          {
            tmp1 = alpha*y[j];
            tmp2 = alpha*x[j];
            for( i = 0; i <= j; ++i )
            { A(i,j) += x[i]*tmp1 + y[i]*tmp2; }
          }
        }
      }
      else
      {
        for( j = 0; j < n; ++j )
        {
          if( (0.0 != x[jx]) || (0.0 != y[jy]) )
          {
            tmp1 = alpha*y[jy];
            tmp2 = alpha*x[jx];
            ix = kx;
            iy = ky;
            for( i = 0; i <= j; ++i )
            {
              A(i,j) += x[ix]*tmp1 + y[iy]*tmp2;
              ix += incx;
              iy += incy;
            }
          }
          jx += incx;
          jy += incy;
        }
      }
    }
    else
    {
      // Form  A  when A is stored in the lower triangle.

      if( (1 == incx) && (1 == incy) )
      {
        for( j = 0; j < n; ++j )
        {
          if( (0.0 != x[j]) || (0.0 != y[j]) )
          {
            tmp1 = alpha*y[j];
            tmp2 = alpha*x[j];
            for( i = j; i < n; ++i )
            { A(i,j) += x[i]*tmp1 + y[i]*tmp2; }
          }
        }
      }
      else
      {
        for( j = 0; j < n; ++j )
        {
          if( (0.0 != x[jx]) || (0.0 != y[jy]) )
          {
            tmp1 = alpha*y[jy];
            tmp2 = alpha*x[jx];
            ix = jx;
            iy = jy;
            for( i = j; i < n; ++i )
            {
              A(i,j) += x[ix]*tmp1 + y[iy]*tmp2;
              ix += incx;
              iy += incy;
            }
          }
          jx += incx;
          jy += incy;
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // dsyr2k
  //------------------------------------------------------------------------

  void dsyr2k( HalfOpt uplo, TransposeOpt trans, int n, int k, double alpha, double *pA, int lda, double *pB, int ldb, double beta, double *pC, int ldc )
  {
    double tmp1, tmp2;
    int i, j, h, nrowa;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };
    auto B = [&]( int i, int j ) noexcept -> double &
    { return pB[col_mjr(i,j,ldb)]; };
    auto C = [&]( int i, int j ) noexcept -> double &
    { return pC[col_mjr(i,j,ldc)]; };

    // Test the input parameters.

    if( kNoTrans == trans )
    { nrowa = n; }else
    { nrowa = k; }

    if( n < 0 ){ throw XerblaException{ "dsyr2k", 2 }; }
    if( k < 0 ){ throw XerblaException{ "dsyr2k", 3 }; }
    if( lda < max(1,nrowa) ){ throw XerblaException{ "dsyr2k", 7 }; }
    if( ldb < max(1,nrowa) ){ throw XerblaException{ "dsyr2k", 9 }; }
    if( ldc < max(1,n) ){ throw XerblaException{ "dsyr2k", 12 }; }

    // Quick return if possible.

    if( (0 == n) || (((0.0 == alpha) || (0 == k)) && (1.0 == beta) )){ return; }

    if( 0.0 == alpha )
    {
      if( kUpper == uplo )
      {
        if( 0.0 == beta )
        {
          for( j = 0; j < n; ++j )
          { for( i = 0; i <= j; ++i )
          { C(i,j) = 0.0; } }
        }
        else
        {
          for( j = 0; j < n; ++j )
          { for( i = 0; i <= j; ++i )
          { C(i,j) = beta*C(i,j); } }
        }
      }
      else
      {
        if( 0.0 == beta )
        {
          for( j = 0; j < n; ++j )
          { for( i = j; i < n; ++i )
          { C(i,j) = 0.0; } }
        }
        else
        {
          for( j = 0; j < n; ++j )
          { for( i = 0; i < n; ++i )
          { C(i,j) += beta*C(i,j); } }
        }
      }
      return;
    }

    // Start the operations.

    if( kNoTrans == trans )
    {
      // Form  C := alpha*A*(~B) + alpha*B*(~A) + C.

      if( kUpper == uplo )
      {
        for( j = 0; j < n; ++j )
        {
          if( 0.0 == beta )
          { for( i = 0; i <= j; ++i )
          { C(i,j) = 0.0; } }
          else if( 1.0 != beta )
          { for( i = 0; i <= j; ++i )
          { C(i,j) = beta*C(i,j); } }
          for( h = 0; h <= j; ++h )
          {
            if( (0.0 != A(j,h)) || (0.0 != B(j,h)) )
            {
              tmp1 = alpha*B(j,h);
              tmp2 = alpha*A(j,h);
              for( i = 0; i <= j; ++i )
              { C(i,j) += A(i,h)*tmp1 + B(i,h)*tmp2; }
            }
          }
        }
      }
      else
      {
        for( j = 0; j < n; ++j )
        {
          if( 0.0 == beta )
          { for( i = j; i < n; ++i )
          { C(i,j) = 0.0; } }
          else if( 1.0 != beta )
          { for( i = j; i < n; ++i )
          { C(i,j) = beta*C(i,j); } }
        }
        for( h = 0; h < k; ++h )
        {
          if( (0.0 != A(j,h)) || (0.0 != B(j,h)) )
          {
            tmp1 = alpha*B(j,h);
            tmp2 = alpha*A(j,h);
            for( i = j; i < n; ++i )
            { C(i,j) += A(i,h)*tmp1 + B(i,h)*tmp2; }
          }
        }
      }
    }
    else
    {
      // Form  C := alpha*(~A)*B + alpha*(~B)*A + C.

      if( kUpper == uplo )
      {
        for( j = 0; j < n; ++j )
        {
          for( i = 0; i <= j; ++i )
          {
            tmp1 = 0.0;
            tmp2 = 0.0;
            for( h = 0; h < k; ++h )
            {
              tmp1 += A(h,i)*B(h,j);
              tmp2 += B(h,i)*A(h,j);
            }
            if( 0.0 == beta )
            { C(i,j) = alpha*tmp1 + alpha*tmp2; }else
            { C(i,j) = beta*C(i,j) + alpha*tmp1 + alpha*tmp2; }
          }
        }
      }
      else
      {
        for( j = 0; j < n; ++j )
        {
          for( i = j; i < n; ++i )
          {
            tmp1 = 0.0;
            tmp2 = 0.0;
            for( h = 0; h < k; ++h )
            {
              tmp1 += A(h,i)*B(h,j);
              tmp2 += B(h,i)*A(h,j);
            }
            if( 0.0 == beta )
            { C(i,j) = alpha*tmp1 + alpha*tmp2; }else
            { C(i,j) = beta*C(i,j) + alpha*tmp1 + alpha*tmp2; }
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // dnrm2
  //------------------------------------------------------------------------

  double dnrm2( int n, double *x, int incx )
  {
    double absxi, norm, scale, ssq;
    int ix;

    if( (n < 1) || (incx < 1) )
    {
      norm = 0.0;
    }
    else if( 1 == n )
    {
      norm = fabs(x[0]);
    }
    else
    {
      scale = 0.0;
      ssq = 1.0;

      // The following loop is equivalent to this call to the LAPACK
      // auxiliary routine:
      // CALL DLASSQ( N, X, INCX, SCALE, SSQ )

      for( ix = 0; ix < (1 + (n-1)*incx); ix += incx )
      {
        if( 0.0 != x[ix] )
        {
          absxi = fabs(x[ix]);
          if( scale < absxi )
          {
            ssq = 1.0 + ssq*sqr(scale/absxi);
            scale = absxi;
          }
          else
          {
            ssq += sqr(absxi/scale);
          }
        }
        norm = scale*sqrt(ssq);
      }
    }

    return norm;
  }

  //------------------------------------------------------------------------
  // dger
  //------------------------------------------------------------------------

  void dger( int m, int n, double alpha, double *x, int incx, double *y, int incy, double *pA, int lda )
  {
    double tmp;
    int i, ix, j, jy, kx;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Test the input parameters.

    if( m < 0 ){ throw XerblaException{ "dger", 1 }; }
    if( n < 0 ){ throw XerblaException{ "dger", 2 }; }
    if( 0 == incx ){ throw XerblaException{ "dger", 5 }; }
    if( 0 == incy ){ throw XerblaException{ "dger", 7 }; }
    if( lda < max( 1, m ) ){ throw XerblaException{ "dger", 9 }; }

    // Quick return if possible.

    if( (0 == m) || (0 == n) || (0.0 == alpha) ){ return; }

    /// Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through A.

    if( incy > 0 )
    { jy = 0; }else
    { jy = (1-(n-1)*incy)-1; }
    if( 1 == incx )
    {
      for( j = 0; j < n; ++j )
      {
        tmp = alpha*y[jy];
        for( i = 0; i < m; ++i )
        { A(i,j) += x[i]*tmp; }
        jy += incy;
      }
    }
    else
    {
      if( incx > 0 )
      { kx = 1; }else
      { kx = (1-(m-1)*incx)-1; }
      for( j = 0; j < n; ++j )
      {
        if( 0.0 != y[jy] )
        {
          tmp = alpha*y[jy];
          ix = kx;
          for( i = 0; i < m; ++i )
          {
            A(i,j) += x[ix]*tmp;
            ix += incx;
          }
        }
        jy += incy;
      }
    }
  }

  //------------------------------------------------------------------------
  // dgemv
  //------------------------------------------------------------------------

  void dgemv( TransposeOpt trans, int m, int n, double alpha, double *pA, int lda, double *x, int incx, double beta, double *y, int incy )
  {
    double tmp;
    int i, ix, iy, j, jx, jy, kx, ky, lenx, leny;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Test the input parameters.

    if( m < 0 ){ throw XerblaException{ "dgemv", 2 }; }
    if( n < 0 ){ throw XerblaException{ "dgemv", 3 }; }
    if( lda < max( 1, m ) ){ throw XerblaException{ "dgemv", 6 }; }
    if( 0 == incx ){ throw XerblaException{ "dgemv", 8 }; }
    if( 0 == incy ){ throw XerblaException{ "dgemv", 11 }; }

    // Quick return if possible.

    if( (0 == m) || (0  == n) ||
      (0.0 == alpha) || (1.0 == beta) ){ return; }

    // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    // up the start points in  X  and  Y.

    if( kNoTrans == trans )
    {
      lenx = n;
      leny = m;
    }
    else
    {
      lenx = m;
      leny = n;
    }
    if( incx > 0 )
    { kx = 0; }else
    { kx = (1-(lenx-1)*incx)-1; }
    if( incy > 0 )
    { ky = 0; }else
    { ky = (1-(leny-1)*incy)-1; }

    // Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through A.
    //
    // First form  y := beta*y.

    if( 1.0 != beta )
    {
      if( 1 == incy )
      {
        if( 0.0 == beta )
        {
          for( i = 0; i < leny; ++i )
          { y[i] = 0.0; }
        }
        else
        {
          for( i = 0; i < leny; ++i )
          { y[i] = beta*y[i]; }
        }
      }
      else
      {
        iy = ky;
        if( 0.0 == beta )
        {
          for( i = 0; i < leny; ++i )
          {
            y[iy] = 0.0;
            iy += incy;
          }
        }
        else
        {
          for( i = 0; i < leny; ++i )
          {
            y[iy] = beta*y[iy];
            iy += incy;
          }
        }
      }
    }

    if( 0.0 == alpha ){ return; }

    if( kNoTrans == trans )
    {
      // Form  y := alpha*A*x + y.

      jx = kx;
      if( 1 == incy )
      {
        for( j = 0; j < n; ++j )
        {
          tmp = alpha*x[jx];
          for( i = 0; i < m; ++i )
          { y[i] += tmp*A(i,j); }
          jx += incx;
        }
      }
      else
      {
        for( j = 0; j < n; ++j )
        {
          tmp = alpha*x[jx];
          iy = ky;
          for( i = 0; i < m; ++i )
          {
            y[iy] += tmp*A(i,j);
            iy += incy;
          }
          jx += incx;
        }
      }
    }
    else
    {
      // Form  y := alpha*(~A)*x + y.

      jy = ky;
      if( 1 == incx )
      {
        for( j = 0; j < n; ++j )
        {
          tmp = 0.0;
          for( i = 0; i < m; ++i )
          { tmp += A(i,j)*x[i]; }
          y[jy] += alpha*tmp;
          jy += incy;
        }
      }
      else
      {
        for( j = 0; j < n; ++j )
        {
          tmp = 0.0;
          ix = kx;
          for( i = 0; i < m; ++i )
          {
            tmp += A(i,j)*x[ix];
            ix += incx;
          }
          y[jy] += alpha*tmp;
          jy += incy;
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // dgemm
  //------------------------------------------------------------------------

  void dgemm( TransposeOpt transa, TransposeOpt transb, int m, int n, int k,
    double alpha, double *pA, int lda, double *pB, int ldb, double beta, double *pC, int ldc )
  {
    double tmp;
    int i, j, h, ncola, nrowa, nrowb;
    bool nota, notb;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };
    auto B = [&]( int i, int j ) noexcept -> double &
    { return pB[col_mjr(i,j,ldb)]; };
    auto C = [&]( int i, int j ) noexcept -> double &
    { return pC[col_mjr(i,j,ldc)]; };

    // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    // transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
    // and  columns of  A  and the  number of  rows  of  B  respectively.

    nota = ( kNoTrans == transa );
    notb = ( kNoTrans == transb );

    if( nota )
    {
      nrowa = m;
      ncola = k;
    }
    else
    {
      nrowa = k;
      ncola = m;
    }
    if( notb )
    { nrowb = k; }else
    { nrowb = n; }

    // Test the input parameters.

    if( m < 0 ){ throw XerblaException{ "dgemm", 3 }; }
    if( n < 0 ){ throw XerblaException{ "dgemm", 4 }; }
    if( k < 0 ){ throw XerblaException{ "dgemm", 5 }; }
    if( lda < max(1,nrowa) ){ throw XerblaException{ "dgemm", 8 }; }
    if( ldb < max(1,nrowb) ){ throw XerblaException{ "dgemm", 10 }; }
    if( ldc < max(1,m) ){ throw XerblaException{ "dgemm", 13 }; }

    // Quick return if possible.

    if( (0 == m) || (0 == n) ||
      (((0.0 == alpha) || (0 == k)) && (1.0 == beta)) ){ return; }

    if( 0.0 == alpha )
    {
      if( 0.0 == beta )
      {
        for( j = 0; j < n; ++j )
        { for( i = 0; i < m; ++i )
        { C(i,j) = 0.0; } }
      }
      else
      {
        for( j = 0; j < n; ++j )
        { for( i = 0; i < m; ++i )
        { C(i,j) = beta*C(i,j); } }
      }
      return;
    }

    // Start the operations.

    if( notb )
    {
      if( nota )
      {
        // Form  C := alpha*A*B + beta*C.

        for( j = 0; j < n; ++j )
        {
          if( 0.0 == beta )
          {
            for( i = 0; i < m; ++i )
            { C(i,j) = 0.0; }
          }
          else if( 1.0 != beta )
          {
            for( i = 0; i < m; ++i )
            { C(i,j) = beta*C(i,j); }
          }
          for( h = 0; h < k; ++h )
          {
            tmp = alpha*B(h,j);
            for( i = 0; i < m; ++i )
            { C(i,j) += tmp*A(i,h); }
          }
        }
      }
      else
      {
        // Form  C := alpha*(~A)*B + beta*C

        for( j = 0; j < n; ++j )
        {
          for( i = 0; i < m; ++i )
          {
            tmp = 0.0;
            for( h = 0; h < k; ++h )
            { tmp += A(h,i)*B(h,j); }
            if( 0.0 == beta )
            { C(i,j) = alpha*tmp; }else
            { C(i,j) = alpha*tmp + beta*C(i,j); }
          }
        }
      }
    }
    else
    {
      if( nota )
      {
        // Form  C := alpha*A*(~B) + beta*C

        for( j = 0; j < n; ++j )
        {
          if( 0.0 == beta )
          {
            for( i = 0; i < m; ++i )
            { C(i,j) = 0.0; }
          }
          else if( 1.0 != beta )
          {
            for( i = 0; i < m; ++i )
            { C(i,j) = beta*C(i,j); }
          }
          for( h = 0; h < k; ++h )
          {
            tmp = alpha*B(j,h);
            for( i = 0; i < m; ++i )
            { C(i,j) += tmp*A(i,h); }
          }
        }
      }
      else
      {
        // Form  C := alpha*(~A)*(~B) + beta*C

        for( j = 0; j < n; ++j )
        {
          for( i = 0; i < m; ++i )
          {
            tmp = 0.0;
            for( h = 0; h < k; ++h )
            { tmp += A(h,i)*B(j,h); }
            if( 0.0 == beta )
            { C(i,j) = alpha*tmp; }else
            { C(i,j) = alpha*tmp + beta*C(i,j); }
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // dlasnst
  //------------------------------------------------------------------------

  double dlanst( NormOpt norm, int n, double *d, double *e )
  {
    int i;
    double anorm, scale, sum;

    if( n <= 0 )
    {
      anorm = 0.0;
    }
    else if( kMaxNorm == norm )
    {
      // Find max(fabs(A(i,j))).

      anorm = fabs(d[n-1]);
      for( i = 0; i < n-1; ++i )
      {
        sum = fabs(d[i]);
        if( (anorm < sum) || isnan(sum) ){ anorm = sum; }
        sum = fabs(e[i]);
        if( (anorm < sum) || isnan(sum) ){ anorm = sum; }
      }
    }
    else if( ( kOneNorm == norm ) || ( kInfNorm == norm ) )
    {
      // Find norm1(A).

      if( 1 == n )
      {
        anorm = fabs(d[0]);
      }
      else
      {
        anorm = fabs(d[0]) + fabs(e[0]);
        sum = fabs(e[n-2])+fabs(d[n-1]);
        if( (anorm < sum) || isnan(sum) ){ anorm = sum; }
        for( i = 1; i < (n-1); ++i )
        {
          sum = fabs(d[i]) + fabs(e[i]) + fabs(e[i-1]);
          if( (anorm < sum) || isnan(sum) ){ anorm = sum; }
        }
      }
    }
    else if( kFrobNorm == norm )
    {
      // Find normF(A).

      scale = 0.0;
      sum = 1.0;
      if( n > 1 )
      {
        dlassq( n-1, e, 1, scale, sum );
        sum = 2*sum;
      }
      dlassq( n, d, 1, scale, sum );
      anorm = scale*sqrt( sum );
    }

    return anorm;
  }

  //------------------------------------------------------------------------
  // dlasnsy
  //------------------------------------------------------------------------

  double dlansy( NormOpt norm, HalfOpt uplo, int n, double *pA, int lda, double *work )
  {
    int i, j;
    double absa, sum, value;
    double ssq[2], colssq[2];

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    if( 0 == n )
    {
      value = 0.0;
    }
    else if( kMaxNorm == norm )
    {
      // Find max(fabs(A(i,j))).

      value = 0.0;
      if( kUpper == uplo )
      {
        for( j = 0; j < n; ++j )
        {
          for( i = 0; i <= j; ++i )
          {
            sum = fabs(A(i,j));
            if( ( value < sum ) || isnan(sum) ){ value = sum; }
          }
        }
      }
      else
      {
        for( j = 0; j < n; ++j )
        {
          for( i = j; i < n; ++i )
          {
            sum = fabs(A(i,j));
            if( ( value < sum ) || isnan(sum) ){ value = sum; }
          }
        }
      }
    }
    else if( (kInfNorm == norm) || (kOneNorm == norm) )
    {
      // Find normI(A) ( = norm1(A), since A is symmetric).

      value = 0.0;
      if( kUpper == uplo )
      {
        for( j = 0; j < n; ++j )
        {
          sum = 0.0;
          for( i = 0; i <= (j-1); ++i )
          {
            absa = fabs(A(i,j));
            sum += absa;
            work[i] += absa;
          }
          work[j] = sum + fabs(A(j,j));
        }
        for( i = 0; i < n; ++i )
        {
          sum = work[i];
          if( (value < sum) || isnan(sum) ){ value = sum; }
        }
      }
      else
      {
        for( i = 0; i < n; ++i )
        { work[i] = 0.0; }
        for( j = 0; j < n; ++j )
        {
          sum = work[j] + fabs(A(j,j));
          for( i = j+1; i < n; ++i )
          {
            absa = fabs(A(i,j));
            sum += absa;
            work[i] += absa;
          }
          if( (value < sum) || isnan(sum) ){ value = sum; }
        }
      }
    }
    else if( kFrobNorm == norm )
    {
      // Find normF(A).
      // SSQ(1) is scale
      // SSQ(2) is sum-of-squares
      // For better accuracy, sum each column separately.

      ssq[0] = 0.0;
      ssq[1] = 1.0;

      // Sum off-diagonals

      if( kUpper == uplo )
      {
        for( j = 1; j < n; ++j )
        {
          colssq[0] = 0.0;
          colssq[1] = 1.0;
          dlassq( (j-1)+1, &A(0,j), 1, colssq[0], colssq[1] );
          dcombssq( ssq, colssq );
        }
      }
      else
      {
        for( j = 0; j < n-1; ++j )
        {
          colssq[0] = 0.0;
          colssq[1] = 1.0;
          dlassq( n-(j+1), &A(j+1,j), 1, colssq[0], colssq[1] );
          dcombssq( ssq, colssq );
        }
      }

      ssq[1] = 2*ssq[1];

      // Sum diagonal

      colssq[0] = 0.0;
      colssq[1] = 1.0;
      dlassq( n, pA, lda+1, colssq[0], colssq[1] );
      dcombssq( ssq, colssq );
      value = ssq[0]*sqrt( ssq[1] );
    }

    return value;
  }


  //------------------------------------------------------------------------
  // dtrmv
  //------------------------------------------------------------------------

  void dtrmv( HalfOpt uplo, TransposeOpt trans, DiagIsUnitOpt diag, 
    int n, double *pA, int lda, double *x, int incx )
  {
    double tmp;
    int i, ix, j, jx, kx;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Test the input parameters.

    if( n < 0 ){ throw XerblaException{ "dtrmv", 4 }; }
    if( lda < max(1,n) ){ throw XerblaException{ "dtrmv", 6 }; }
    if( 0 == incx ){ throw XerblaException{ "dtrmv", 8 }; }

    // Quick return if possible.
    if( 0 == n ){ return; }

    if( incx < 0 )
    { kx = (1-(n-1)*incx)-1; }else
    { kx = 0; }

    // Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through A.

    if( kNoTrans == trans )
    {
      // Form  x := A*x.

      if( kUpper == uplo )
      {
        if( 1 == incx )
        {
          for( j = 0; j < n; ++j )
          {
            if( 0.0 != x[j] )
            {
              tmp = x[j];
              for( i = 0; i <= j-1; ++i )
              { x[i] += tmp*A(i,j); }
              if( kDiagNotUnit == diag )
              { x[j] *= A(j,j); }
            }
          }
        }
        else
        {
          jx = kx;
          for( j = 0; j < n; ++j )
          {
            if( 0.0 != x[j] )
            {
              tmp = x[j];
              ix = kx;
              for( i = 0; i <= j-1; ++i )
              {
                x[ix] += tmp*A(i,j);
                ix += incx;
              }
              if( kDiagNotUnit == diag )
              { x[j] *= A(j,j); }
            }
            jx += incx;
          }
        }
      }
      else
      {
        if( 1 == incx )
        {
          for( j = n-1; j >= 0; --j )
          {
            if( 0.0 != x[j] )
            {
              tmp = x[j];
              for( i = n-1; i >= j+1; --i )
              { x[i] += tmp*A(i,j); }
              if( kDiagNotUnit == diag )
              { x[j] *= A(j,j); }
            }
          }
        }
        else
        {
          kx += (n-1)*incx;
          jx = kx;
          for( j = n-1; j >= 0; --j )
          {
            if( 0.0 != x[jx] )
            {
              tmp = x[jx];
              ix = kx;
              for( i = n-1; i >= j+1; --i )
              {
                x[ix] += tmp*A(i,j);
                ix += incx;
              }
              if( kDiagNotUnit == diag )
              { x[jx] *= A(j,j); }
            }
            jx += incx;
          }
        }
      }
    }
    else
    {
      // Form  x := (~A)*x.

      if( kUpper == uplo )
      {
        if( 1 == incx )
        {
          for( j = n-1; j >= 0; --j )
          {
            tmp = x[j];
            if( kDiagNotUnit == diag )
            { tmp *= A(j,j); }
            for( i = j-1; i >= 0; --i )
            { tmp += A(i,j)*x[i]; }
            x[j] = tmp;
          }
        }
        else
        {
          jx = kx + (n-1)*incx;

          for( j = n-1; j >= 0; --j )
          {
            tmp = x[jx];
            ix = jx;
            if( kDiagNotUnit == diag )
            { tmp *= A(j,j); }
            for( i = j-1; i >= 0; --i )
            {
              ix -= incx;
              tmp += A(i,j)*x[ix];
            }
            x[jx] = tmp;
            jx -= incx;
          }
        }
      }
      else
      {
        if( 1 == incx )
        {
          for( j = 0; j < n; ++j )
          {
            tmp = x[j];
            if( kDiagNotUnit == diag )
            { tmp *= A(j,j); }
            for( i = j+1; i < n; ++i )
            { tmp += A(i,j)*x[i]; }
            x[j] = tmp;
          }
        }
        else
        {
          jx = kx;
          for( j = 0; j < n; ++j )
          {
            tmp = x[jx];
            ix = jx;
            if( kDiagNotUnit == diag )
            { tmp *= A(j,j); }
            for( i = j+1; i < n; ++i )
            {
              ix += incx;
              tmp += A(i,j)*x[ix];
            }
            x[jx] = tmp;
            jx += incx;
          }
        }
      }
    }
  }


  //------------------------------------------------------------------------
  // dtrmm
  //------------------------------------------------------------------------

  void dtrmm( SideOpt side, HalfOpt uplo, TransposeOpt transa, DiagIsUnitOpt diag,
    int m, int n, double alpha, double *pA, int lda, double *pB, int ldb )
  {
    double tmp;
    int i, j, k, nrowa;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };
    auto B = [&]( int i, int j ) noexcept -> double &
    { return pB[col_mjr(i,j,ldb)]; };

    // Test the input parameters.

    if( kLeft == side )
    { nrowa = m; }else
    { nrowa = n; }

    if( m < 0 ){ throw XerblaException{ "dtrmm", 5 }; }
    if( n < 0 ){ throw XerblaException{ "dtrmm", 6 }; }
    if( lda < max(1,nrowa) ){ throw XerblaException{ "dtrmm", 9 }; }
    if( ldb < max(1,m) ){ throw XerblaException{ "dtrmm", 11 }; }

    // Quick return if possible.

    if( (0 == m) || (0 == n) ){ return; }

    if( kLeft == side )
    {
      if( kNoTrans == transa )
      {
        // Form  B := alpha*A*B.

        if( kUpper == uplo )
        {
          for( j = 0; j < n; ++j )
          {
            for( k = 0; k < m; ++k )
            {
              if( 0.0 != B(k,j) )
              {
                tmp = alpha*B(k,j);
                for( i = 0; i <= k-1; ++i )
                { B(i,j) += tmp*A(i,k); }
                if( kDiagNotUnit == diag )
                { tmp = tmp*A(k,k); }
                B(k,j) = tmp;
              }
            }
          }
        }
        else
        {
          for( j = 0; j < n; ++j )
          {
            for( k = m-1; k >= 0; --k )
            {
              if( 0.0 != B(k,j) )
              {
                tmp = alpha*B(k,j);
                B(k,j) = tmp;
                if( kDiagNotUnit == diag )
                { B(k,j) *= A(k,k); }
                for( i = k+1; i < m; ++i )
                { B(i,j) += tmp*A(i,k); }
              }
            }
          }
        }
      }
      else
      {
        // Form  B := alpha*(~A)*B.

        if( kUpper == uplo )
        {
          for( j = 0; j < n; ++j )
          {
            for( i = m-1; i >= 0; --i )
            {
              tmp = B(i,j);
              if( kDiagNotUnit == diag )
              { tmp *= A(i,i); }
              for( k = 0; k <= i-1; ++k )
              { tmp += A(k,i)*B(k,j); }
              B(i,j) = alpha*tmp;
            }
          }
        }
        else
        {
          for( j = 0; j < n; ++j )
          {
            for( i = 0; i < m; ++i )
            {
              tmp = B(i,j);
              if( kDiagNotUnit == diag )
              { tmp *= A(i,i); }
              for( k = i+1; k < m; ++k )
              { tmp += A(k,i)*B(k,j); }
              B(i,j) = alpha*tmp;
            }
          }
        }
      }
    }
    else
    {
      if( kNoTrans == transa )
      {
        // Form  B := alpha*B*A.

        if( kUpper == uplo )
        {
          for( j = n-1; j >= 0; --j )
          {
            tmp = alpha;
            if( kDiagNotUnit == diag )
            { tmp *= A(j,j); }
            for( i = 0; i < m; ++i )
            { B(i,j) = tmp*B(i,j); }
            for( k = 0; k <= j-1; ++k )
            {
              if( 0.0 != A(k,j) )
              {
                tmp = alpha*A(k,j);
                for( i = 0; i < m; ++i )
                { B(i,j) += tmp*B(i,k); }
              }
            }
          }
        }
        else
        {
          for( j = 0; j < n; ++j )
          {
            tmp = alpha;
            if( kDiagNotUnit == diag )
            { tmp *= A(j,j); }
            for( i = 0; i < m; ++i )
            { B(i,j) = tmp*B(i,j); }
            for( k = j+1; k < n; ++k )
            {
              if( 0.0 != A(k,j) )
              {
                tmp = alpha*A(k,j);
                for( i = 0; i < m; ++i )
                { B(i,j) += tmp*B(i,k); }
              }
            }
          }
        }
      }
      else
      {
        // Form  B := alpha*B*(~A).

        if( kUpper == uplo )
        {
          for( k = 0; k < n; ++k )
          {
            for( j = 0; j <= k-1; ++j )
            {
              if( 0.0 != A(j,k) )
              {
                tmp = alpha*A(j,k);
                for( i = 0; i < m; ++i )
                { B(i,j) += tmp*B(i,k); }
              }
            }
            tmp = alpha;
            if( kDiagNotUnit == diag )
            { tmp *= A(k,k); }
            if( 1.0 != tmp )
            {
              for( i = 0; i < m; ++i )
              { B(i,k) = tmp*B(i,k); }
            }
          }
        }
        else
        {
          for( k = n-1; k >= 0; --k )
          {
            for( j = k+1; j < n; ++j )
            {
              if( 0.0 != A(j,k) )
              {
                tmp = alpha*A(j,k);
                for( i = 0; i < m; ++i )
                { B(i,j) += tmp*B(i,k); }
              }
            }
            tmp = alpha;
            if( kDiagNotUnit == diag )
            { tmp *= A(k,k); }
            if( 1.0 != tmp )
            {
              for( i = 0; i < m; ++i )
              { B(i,k) = tmp*B(i,k); }
            }
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // dlarf
  //------------------------------------------------------------------------

  void dlarf( SideOpt side, int m, int n, double *v, int incv, double tau, double *pC, int ldc, double *work )
  {
    bool applyleft;
    int i, lastv, lastc;

    applyleft = ( kLeft == side );
    lastv = 0;
    lastc = 0;

    // Set up variables for scanning V.  LASTV begins pointing to the end
    // of V.
    if( 0.0 != tau )
    {
      if( applyleft )
      { lastv = m; }else
      { lastv = n; }
      if( incv > 0 )
      { i = (1+(lastv-1)*incv)-1; }else
      { i = 0; }

      // Look for the last non-zero row in V.
      while( (lastv > 0) && (0.0 == v[i]) )
      { --lastv; i -= incv; }

      // Scan for the last non-zero column in C(0:lastv-1,:).
      if( applyleft )
      { lastc = iladlc( lastv, n, pC, ldc ) + 1; }
      // Scan for the last non-zero row in C(:,0:lastv-1).
      else
      { lastc = iladlr( m, lastv, pC, ldc ) + 1; }

      // Note that lastc == 0 renders the BLAS operations null; no special
      // case is needed at this level.
      if( applyleft )
      {
        // Form H*C
        if( lastv > 0 )
        {
          dgemv( kTrans, lastv, lastc, 1.0, pC, ldc, v, incv, 0.0, work, 1 );
          dger( lastv, lastc, -tau, v, incv, work, 1, pC, ldc );
        }
      }
      else
      {
        // Form C*H
        if( lastv > 0 )
        {
          dgemv( kNoTrans, lastc, lastv, 1.0, pC, ldc, v, incv, 0.0, work, 1 );
          dger( lastc, lastv, -tau, work, 1, v, incv, pC, ldc );
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // dlarfb
  //------------------------------------------------------------------------

  void dlarfb( SideOpt side, TransposeOpt trans, DirectOpt direct, StoreOpt storev, int m, int n, int k,
    double *pV, int ldv, double *pT, int ldt, double *pC, int ldc, double *pWork, int ldwork )
  {
    TransposeOpt transt;
    int i, j;

    auto V = [&]( int i, int j ) noexcept -> double &
    { return pV[col_mjr(i,j,ldv)]; };
    auto T = [&]( int i, int j ) noexcept -> double &
    { return pT[col_mjr(i,j,ldt)]; };
    auto C = [&]( int i, int j ) noexcept -> double &
    { return pC[col_mjr(i,j,ldc)]; };
    auto Work = [&]( int i, int j ) noexcept -> double &
    { return pWork[col_mjr(i,j,ldwork)]; };

    // Quick return if possible

    if( (m <= 0) || (n <= 0) ){ return; }

    transt = ( kNoTrans == trans ) ? kTrans : kNoTrans;

    if( kByCol == storev )
    {
      if( kForwd == direct )
      {
        // Let  V = [ V1 ]    (first K rows)
        //          [ V2 ]
        //
        //  where V1 is unit lower triangular.

        if( kLeft == side )
        {
          // Form  H*C  or  (~H)*C  where  C = [ C1 ]
          //                                   [ C2 ]
          //
          // W := (~C)*V = ((~C1)*V1 + (~C2)*V2)  (stored in WORK)
          //
          // W := (~C1)

          for( j = 0; j < k; ++k )
          { dcopy( n, &C(j,0), ldc, &Work(0,j), 1 ); }

          // W := W*V1

          dtrmm( kRight, kLower, kNoTrans, kDiagUnit, n, k, 1.0, pV, ldv, pWork, ldwork );

          if( m > k )
          {
            // W := W + (~C2)*V2
            dgemm( kTrans, kNoTrans, n, k, m-k, 1.0, &C((k+1)-1,0), ldc, &V((k+1)-1,0), ldv, 1.0, pWork, ldwork );

            // W := W*(~T)  or  W*T
            dtrmm( kRight, kUpper, transt, kDiagNotUnit, n, k, 1.0, pT, ldt, pWork, ldwork );

            // C := C - V*(~W)
            if( m > k )
            { dgemm( kNoTrans, kTrans, m-k, n, k, -1.0, &V((k+1)-1,0), ldv, pWork, ldwork, 1.0, &C((k+1)-1,0), ldc ); }

            // W := W*(~V1)
            dtrmm( kRight, kLower, kTrans, kDiagUnit, n, k, 1.0, pV, ldv, pWork, ldwork );

            // C1 := C1 - (~W)
            for( j = 0; j < k; ++j )
            { for( i = 0; i < n; ++i )
            { C(j,i) -= Work(i,j); } }
          }
        }
        else if( kRight == side )
        {
          // Form  C*H  or  C*(~H) where  C = [ C1  C2 ]
          //
          // W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
          //
          // W := C1

          for( j = 0; j < k; ++j )
          { dcopy( m, &C(0,j), 1, &Work(0,j), 1 ); }

          // W := W*V1
          dtrmm( kRight, kLower, kNoTrans, kDiagUnit, m, k, 1.0, pV, ldv, pWork, ldwork );

          // W := W + C2*V2
          if( n > k )
          { dgemm( kNoTrans, kNoTrans, m, k, n-k, 1.0, &C(0,(k+1)-1), ldc, &V((k+1)-1,0), ldv, 1.0, pWork, ldwork ); }

          // W := W*T  or  W*(~T)
          dtrmm( kRight, kUpper, trans, kDiagNotUnit, m, k, 1.0, pT, ldt, pWork, ldwork );

          // C := C - W*(~V)
          if( n > k )
          { dgemm( kNoTrans, kTrans, m, n-k, k, -1.0, pWork, ldwork, &V((k+1)-1,0), ldv, 1.0, &C(0,(k+1)-1), ldc ); }

          // W := W*(~V1)
          dtrmm( kRight, kLower, kTrans, kDiagUnit, m, k, 1.0, pV, ldv, pWork, ldwork );

          // C1 := C1 - W
          for( j = 0; j < k; ++j )
          { for( i = 0; i < m; ++i )
          { C(i,j) -= Work(i,j); } }
        }
      }
      else
      {
        // Let  V = [ V1 ]
        //          [ V2 ]    (last K rows)
        //
        // where V2 is unit upper triangular.

        if( kLeft == side )
        {
          // Form  H*C  or  (~H)*C  where  C = [ C1 ]
          //                                   [ C2 ]
          //
          // W := (~C)*V = ((~C1)*V1 + (~C2)*V2)  (stored in WORK)

          // W := (~C2)
          for( j = 0; j < k; ++j )
          { dcopy( n, &C((m-k+(j+1))-1,0), ldc, &Work(0,j), 1 ); }

          // W := W*V2
          dtrmm( kRight, kUpper, kNoTrans, kDiagUnit, n, k, 1.0, &V((m-k+1)-1,0), ldv, pWork, ldwork );

          // W := W + (~C1)*V1
          if( m > k )
          { dgemm( kTrans, kNoTrans, n, k, m-k, 1.0, pC, ldc, pV, ldv, 1.0, pWork, ldwork ); }

          // W := W*(~T) or W*T
          dtrmm( kRight, kLower, transt, kDiagNotUnit, n, k, 1.0, pT, ldt, pWork, ldwork );

          // C := C - V*(~W)
          if( m > k )
          { dgemm( kNoTrans, kTrans, m-k, n, k, -1.0, pV, ldv, pWork, ldwork, 1.0, pC, ldc ); }

          // W := W*(~V2)
          dtrmm( kRight, kUpper, kTrans, kDiagUnit, n, k, 1.0, &V((m-k+1)-1,0), ldv, pWork, ldwork );

          // C2 := C2 - (~W)
          for( j = 0; j < k; ++j )
          { for( i = 0; i < n; ++i )
          { C((m-k+(j+1))-1,i) -= Work(i,j); } }
        }
        else if( kRight == side )
        {
          // Form C*H or C*(~H) where C = [ C1  C2 ]
          //
          // W := C*V = (C1*V1 + C2*V2) (stored in WORK)
          //
          // W := C2
          for( j = 0; j < k; ++j )
          { dcopy( m, &C(0,(n-k+(j+1))-1), 1, &Work(0,j), 1 ); }

          // W := W*V2
          dtrmm( kRight, kUpper, kNoTrans, kDiagUnit, m, k, 1.0, &V((n-k+1)-1,0), ldv, pWork, ldwork );

          // W := W + C1*V1
          if( n > k )
          { dgemm( kNoTrans, kNoTrans, m, k, n-k, 1.0, pC, ldc, pV, ldv, 1.0, pWork, ldwork ); }

          // W := W*T or W*(~T)
          dtrmm( kRight, kLower, trans, kDiagNotUnit, m, k, 1.0, pT, ldt, pWork, ldwork );

          // C := C - W*(~V)
          if( n > k )
          { dgemm( kNoTrans, kTrans, m, n-k, k, -1.0, pWork, ldwork, pV, ldv, 1.0, pC, ldc ); }

          // W := W*(~V2)
          dtrmm( kRight, kUpper, kTrans, kDiagUnit, m, k, 1.0, &V((n-k+1)-1,0), ldv, pWork, ldwork );

          // C2 := C2 - W

          for( j = 0; j < k; ++j )
          { for( i = 0; i < m; ++i )
          { C(i,(n-k+(j+1))-1) -= Work(i,j); } }
        }
      }
    }
    else if( kByRow == storev )
    {
      if( kForwd == direct )
      {
        // Let  V = [ V1  V2 ]    (V1: first K columns)
        // where  V1  is unit upper triangular.

        if( kLeft == side )
        {
          // Form  H * C  or  H**T * C  where  C = [ C1 ]
          //                                       [ C2 ]
          //
          // W := (~C)*(~V) = ((~C1)*(~V1) + (~C2)*(~V2)) (stored in WORK)

          // W := (~C1)
          for( j = 0; j < k; ++j )
          { dcopy( n, &C(j,0), ldc, &Work(0,j), 1 ); }

          // W := W*(~V1)
          dtrmm( kRight, kUpper, kTrans, kDiagUnit, n, k, 1.0, pV, ldv, pWork, ldwork );

          // W := W + (~C2)*(~V2)
          if( m > k )
          { dgemm( kTrans, kTrans, n, k, m-k, 1.0, &C((k+1)-1,0), ldc, &V(0,(k+1)-1), ldv, 1.0, pWork, ldwork ); }

          // W := W*(~T) or W*T
          dtrmm( kRight, kUpper, transt, kDiagNotUnit, n, k, 1.0, pT, ldt, pWork, ldwork );

          // C := C - (~V)*(~W)
          if( m > k )
          { dgemm( kTrans, kTrans, m-k, n, k, -1.0, &V(0,(k+1)-1), ldv, pWork, ldwork, 1.0, &C((k+1)-1,0), ldc ); }

          // W := W*V1
          dtrmm( kRight, kUpper, kNoTrans, kDiagUnit, n, k, 1.0, pV, ldv, pWork, ldwork );

          // C1 := C1 - (~W)
          for( j = 1; j < k; ++j )
          { for( i = 0; i < n; ++i )
          { C(j,i) -= Work(i,j); } }
        }
        else if( kRight == side )
        {
          // Form C*H or C*(~H)  where  C = [ C1  C2 ]
          //
          // W := C*(~V) = (C1*(~V1) + C2*(~V2)) (stored in WORK)

          // W := C1
          for( j = 0; j < k; ++j )
          { dcopy( m, &C(0,j), 1, &Work(0,j), 1 ); }

          // W := W*(~V1)
          dtrmm( kRight, kUpper, kTrans, kDiagUnit, m, k, 1.0, pV, ldv, pWork, ldwork );

          // W := W + C2*(~V2)
          if( n > k )
          { dgemm( kNoTrans, kTrans, m, k, n-k, 1.0, &C(0,(k+1)-1), ldc, &V(0,(k+1)-1), ldv, 1.0, pWork, ldwork ); }

          // W := W*T or W*(~T)
          dtrmm( kRight, kUpper, trans, kDiagNotUnit, m, k, 1.0, pT, ldt, pWork, ldwork );

          // C := C - W*V2
          if( n > k )
          { dgemm( kNoTrans, kNoTrans, m, n-k, k, -1.0, pWork, ldwork, &V(0,(k+1)-1), ldv, 1.0, &C(0,(k+1)-1), ldc ); }

          // W := W*V1
          dtrmm( kRight, kUpper, kNoTrans, kDiagUnit, m, k, 1.0, pV, ldv, pWork, ldwork );

          // C1 := C1 - W
          for( j = 0; j < k; ++j )
          { for( i = 0; i < m; ++i )
          { C(i,j) -= Work(i,j); } }
        }
      }
      else
      {
        // Let  V =  [ V1  V2 ]    (V2: last K columns)
        //
        // where  V2  is unit lower triangular.

        if( kLeft == side )
        {
          // Form H*C or (~H)*C where C = [ C1 ]
          //                              [ C2 ]
          //
          // W := (~C)*(~V) = ((~C1)*(~V1) + (~C2)*(~V2)) (stored in WORK)

          // W := (~C2)
          for( j = 0; j < k; ++j )
          { dcopy( n, &C((m-k+(j+1))-1,0), ldc, &Work(0,j), 1 ); }

          // W := W*(~V2)
          dtrmm( kRight, kLower, kTrans, kDiagUnit, n, k, 1.0, &V(0,(m-k+1)-1), ldv, pWork, ldwork );

          // W := W + (~C1)*(~V1)
          if( m > k )
          { dgemm( kTrans, kTrans, n, k, m-k, 1.0, pC, ldc, pV, ldv, 1.0, pWork, ldwork ); }

          // W := W*(~T) or W*T
          dtrmm( kRight, kLower, transt, kDiagNotUnit, n, k, 1.0, pT, ldt, pWork, ldwork );

          // C := C - (~V)*(~W)
          if( m > k )
          { dgemm( kTrans, kTrans, m-k, n, k, -1.0, pV, ldv, pWork, ldwork, 1.0, pC, ldc ); }

          // W := W*V2
          dtrmm( kRight, kLower, kNoTrans, kDiagUnit, n, k, 1.0, &V(0,(m-k+1)-1), ldv, pWork, ldwork );

          // C2 := C2 - (~W)
          for( j = 0; j < k; ++j )
          { for( i = 0; i < n; ++i )
          { C((m-k+(j+1))-1,0) -= Work(i,j); } }
        }
        else if( kRight == side )
        {
          // Form C*H or C*(~H) where C = [ C1  C2 ]
          //
          // W := C*(~V) = (C1*(~V1) + C2*(~V2))  (stored in WORK)

          // W := C2
          for( j = 0; j < k; ++j )
          { dcopy( m, &C(0,(n-k+(j+1))-1), 0, &Work(0,j), 1 ); }

          // W := W*(~V2)
          dtrmm( kRight, kLower, kTrans, kDiagUnit, m, k, 1.0, &V(0,(n-k+1)-1), ldv, pWork, ldwork );

          // W := W + C1*(~V1)
          if( n > k )
          { dgemm( kNoTrans, kTrans, m, k, n-k, 1.0, pC, ldc, pV, ldv, 1.0, pWork, ldwork ); }

          // W := W*T or W*(~T)
          dtrmm( kRight, kLower, trans, kDiagNotUnit, m, k, 1.0, pT, ldt, pWork, ldwork );

          // C := C - W*V1
          if( n > k )
          { dgemm( kNoTrans, kNoTrans, m, n-k, k, -1.0, pWork, ldwork, pV, ldv, 1.0, pC, ldc ); }

          // W := W*V2
          dtrmm( kRight, kLower, kNoTrans, kDiagUnit, m, n, 1.0, &V(0,(n-k+1)-1), ldv, pWork, ldwork );

          // C1 := C1 - W
          for( j = 0; j < k; ++j )
          { for( i = 0; i < m; ++i )
          { C(i,(n-k+(j+1))-1) -= Work(i,j); } }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // dlarfg
  //------------------------------------------------------------------------

  void dlarfg( int n, double &alpha, double *x, int incx, double &tau )
  {
    int j, knt;
    double beta, rsafmn, safmin, xnorm;

    if( n <= 1 )
    {
      tau = 0.0;
      return;
    }

    xnorm = dnrm2( n-1, x, incx );

    if( 0.0 == xnorm )
    {
      // H  =  I
      tau = 0.0;
    }
    else
    {
      // general case

      beta = -copysign( dlapy2( alpha, xnorm ), alpha );
      safmin = DBL_MIN;
      knt = 0;
      if( fabs(beta) < safmin )
      {
        // XNORM, BETA may be inaccurate; scale X and recompute them

        rsafmn = 1.0/safmin;
        do
        {
          ++knt;
          dscal( n-1, rsafmn, x, incx );
          beta *= rsafmn;
          alpha *= rsafmn;
        }
        while( (fabs(beta) < safmin) && (knt < 20) );

        // New BETA is at most 1, at least SAFMIN

        xnorm = dnrm2( n-1, x, incx );
        beta = -copysign( dlapy2( alpha, xnorm ), alpha );
      }
      tau = ( beta - alpha )/beta;
      dscal( n-1, 1.0/(alpha-beta), x, incx );

      // If ALPHA is subnormal, it may lose relative accuracy

      for( j = 0; j < knt; ++j )
      { beta *= safmin; }
      alpha = beta;
    }
  }

  //------------------------------------------------------------------------
  // dlarft
  //------------------------------------------------------------------------

  void dlarft( DirectOpt direct, StoreOpt storev, int n, int k, double *pV, int ldv, double *tau, double *pT, int ldt )
  {
    int i, j, prevlastv, lastv;

    auto V = [&]( int i, int j ) noexcept -> double &
    { return pV[col_mjr(i,j,ldv)]; };
    auto T = [&]( int i, int j ) noexcept -> double &
    { return pT[col_mjr(i,j,ldt)]; };

    // Quick return if possible
    if( 0 == n ){ return; }

    if( kForwd == direct )
    {
      prevlastv = n-1;
      for( i = 0; i < k; ++i )
      {
        prevlastv = max( i, prevlastv );

        if( 0.0 == tau[i] )
        {
          // H(i) = I
          for( j = 0; j <= i; ++j )
          { T(j,i) = 0.0; }
        }
        else
        {
          // general case
          if( kByCol == storev )
          {
            // Skip any trailing zeros.
            for( lastv = n-1; lastv >= i+1; --lastv )
            { if( 0.0 != V(lastv,i) ){ break; } }

            for( j = 0; j <= i-1; ++j )
            { T(j,i) = -tau[i]*V(i,j); }

            j = min( lastv, prevlastv );

            // T(0:i-1,i) := - tau(i) * (~V(i:j,0:i-1)) * V(i:j,i)
            dgemv( kTrans, j-i, (i-1)+1, -tau[i], &V(i+1,0), ldv, &V(i+1,i), 1, 1.0, &T(0,i), 1 );
          }
          else
          {
            // Skip any trailing zeros.
            for( lastv = n-1; lastv >= i+1; --lastv )
            { if( 0.0 != V(i,lastv) ){ break; } }

            for( j = 0; j <= i-1; ++j )
            { T(j,i) = -tau[i]*V(j,i); }

            j = min( lastv, prevlastv );

            // T(0:i-1,i) := - tau(i) * V(0:i-1,i:j) * (~V(i,i:j))
            dgemv( kNoTrans, (i-1)+1, j-i, -tau[i], &V(0,i+1), ldv, &V(i,i+1), ldv, 1.0, &T(0,i), 1 );
          }

          // T(0:i-1,i) := T(0:i-1,1:i-1) * T(0:i-1,i)
          dtrmv( kUpper, kNoTrans, kDiagNotUnit, (i-1)+1, pT, ldt, &T(0,i), 1 );

          T(i,i) = tau[i];
          if( i > 0 )
          { prevlastv = max( prevlastv, lastv ); }else
          { prevlastv = lastv; }
        }
      }
    }
    else
    {
      prevlastv = 0;
      for( i = k-1; i >= 0; --i )
      {
        if( 0.0 == tau[i] )
        {
          // H(i) = I
          for( j = i; j < k; ++j )
          { T(j,i) = 0.0; }
        }
        else
        {
          // general case

          if( i < k )
          {
            if( 'C' == storev )
            {
              // Skip any leading zeros.
              for( lastv = 0; lastv <= i-1; ++lastv )
              { if( 0.0 != V(lastv,i) ){ break; } }

              for( j = i+1; j < k; ++j )
              { T(j,i) = -tau[i]*V((n-k+(i+1))-1,j); }

              j = max( lastv, prevlastv );

              // T(i+1:k-1,i) = -tau(i) * (~V(j:n-k+i-1,i+1:k-1))*V(j:n-k+i-1,i)
              dgemv( kTrans, (n-k+(i+1)-(j+1))-1, k-(i+1), -tau[i], &V(j,i+1), ldv, &V(j,i), 1, 1.0, &T(i+i,i), 1 );
            }
            else
            {
              // Skip any leading zeros.

              for( lastv = 0; lastv <= i-1; ++lastv )
              { if( 0.0 != V(i,lastv) ){ break; } }

              for( j = i+1; j < k; ++j )
              { T(j,i) = -tau[i]*V(j,(n-k+(i+1))-1); }

              j = max( lastv, prevlastv );

              // T(i+1:k-1,i) = -tau(i) * V(i+1:k-1,j:(n-k+(i+1))-1) * (~V(i,j:(n-k+(i+1))-1))
              dgemv( kNoTrans, k-(i+1), (n-k+(i+1)-(j+1))-1, -tau[i], &V(i+1,j), ldv, &V(i,j), ldv, 1.0, &T(i+1,i), 1 );
            }

            // T(i+1:k-1,i) := T(i+1:k-1,i+1:k-1) * T(i+1:k-1,i)
            dtrmv( kLower, kNoTrans, kDiagNotUnit, k-(i+1), &T(i+1,i+1), ldt, &T(i+1,i), 1 );

            if( i > 1 )
            { prevlastv = min( prevlastv, lastv ); }else
            { prevlastv = lastv; }
          }
          T(i,i) = tau[i];
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // dsymv
  //------------------------------------------------------------------------

  void dsymv( HalfOpt uplo, int n, double alpha, double *pA, int lda, double *x, int incx, double beta, double *y, int incy )
  {
    double tmp1, tmp2;
    int i, ix, iy, j, jx, jy, kx, ky;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Test the input parameters.

    if( n < 0 ){ throw XerblaException{ "dsymv", 2 }; }
    if( lda < max(1,n) ){ throw XerblaException{ "dsymv", 5 }; }
    if( 0 == incx ){ throw XerblaException{ "dsymv", 7 }; }
    if( 0 == incy ){ throw XerblaException{ "dsymv", 10 }; }

    // Quick return if possible.

    if( (0 == n) || ((0.0 == alpha) && (1.0 == beta)) ){ return; }

    // Set up the start points in  X  and  Y.

    if( incx > 0 )
    { kx = 0; }else
    { kx = (1-(n-1)*incx)-1; }
    if( incy > 0 )
    { ky = 0; }else
    { ky = (1-(n-1)*incy)-1; }

    // Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through the triangular part
    // of A.
    //
    // First form  y := beta*y.

    if( 1.0 != beta )
    {
      if( 1 == incy )
      {
        if( 0.0 == beta )
        {
          for( i = 0; i < n; ++i )
          { y[i] = 0.0; }
        }
        else
        {
          for( i = 0; i < n; ++i )
          { y[i] = beta*y[i]; }
        }
      }
      else
      {
        iy = ky;
        if( 0.0 == beta )
        {
          for( i = 0; i < n; ++i )
          {
            y[iy] = 0.0;
            iy += incy;
          }
        }
        else
        {
          for( i = 0; i < n; ++i )
          {
            y[iy] = beta*y[iy];
            iy += incy;
          }
        }
      }
    }

    if( 0.0 == alpha ){ return; }

    if( kUpper == uplo )
    {
      // Form  y  when A is stored in upper triangle.

      if( (1 == incx) && (1 == incy) )
      {
        for( j = 0; j < n; ++j )
        {
          tmp1 = alpha*x[j];
          tmp2 = 0.0;
          for( i = 0; i <= (j-1); ++i )
          {
            y[i] += tmp1*A(i,j);
            tmp2 += A(i,j)*x[i];
          }
          y[j] += tmp1*A(j,j) + alpha*tmp2;
        }
      }
      else
      {
        jx = kx;
        jy = ky;
        for( j = 0; j < n; ++j )
        {
          tmp1 = alpha*x[jx];
          tmp2 = 0.0;
          ix = kx;
          iy = ky;
          for( i = 0; i <= (j-1); ++i )
          {
            y[iy] += tmp1*A(i,j);
            tmp2 += A(i,j)*x[ix];
            ix += incx;
            iy += incy;
          }
          y[jy] += tmp1*A(j,j) + alpha*tmp2;
          jx += incx;
          jy += incy;
        }
      }
    }
    else
    {
      // Form  y  when A is stored in lower triangle.

      if( (1 == incx) && (1 == incy) )
      {
        for( j = 0; j < n; ++j )
        {
          tmp1 = alpha*x[j];
          tmp2 = 0.0;
          y[j] += tmp1*A(j,j);
          for( i = j; i < n; ++i )
          {
            y[i] += tmp1*A(i,j);
            tmp2 += A(i,j)*x[i];
          }
          y[j] += alpha*tmp2;
        }
      }
      else
      {
        jx = kx;
        jy = ky;
        for( j = 0; j < n; ++j )
        {
          tmp1 = alpha*x[jx];
          tmp2 = 0.0;
          y[jy] += tmp1*A(j,j);
          ix = jx;
          iy = jy;
          for( i = j+1; i < n; ++i )
          {
            ix += incx;
            iy += incy;
            y[iy] += tmp1*A(i,j);
            tmp2 += A(i,j)*x[ix];
          }
          y[jy] += alpha*tmp2;
          jx += incx;
          jy += incy;
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // dsytd2
  //------------------------------------------------------------------------

  void dsytd2( HalfOpt uplo, int n, double *pA, int lda, double *d, double *e, double *tau )
  {
    int i;
    double alpha, taui;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // Test the input parameters

    if( n < 0 ){ throw XerblaException{ "dsytd2", 2 }; }
    if( lda < max(1,n) ){ throw XerblaException{ "dsytd2", 4 }; }

    // Quick return if possible

    if( n <= 0 ){ return; }

    if( kUpper == uplo )
    {
      // Reduce the upper triangle of A

      for( i = ((n-1)-1); i >= 0; --i )
      {
        // Generate elementary reflector H(i) = I - tau*v*(~v)
        // to annihilate A(0:i-1,i+1)

        dlarfg( i+1, A(i,i+1), &A(0,i+1), 1, taui );
        e[i] = A(i,i+1);

        if( 0.0 != taui )
        {
          // Apply H(i) from both sides to A(0:i,0:i)

          A(i,i+1) = 1.0;

          // Compute  x := tau*A*v  storing x in TAU(0:i)
          dsymv( kUpper, i+1, taui, pA, lda, &A(0,i+1), 1, 0.0, tau, 1 );

          // Compute  w := x - ( (1/2)*tau*dot(x,v) )*v
          alpha = -0.5*taui*ddot( i+1, tau, 1, &A(0,i+1), 1 );
          daxpy( i+1, alpha, &A(0,i+1), 1, tau, 1 );

          // Apply the transformation as a rank-2 update:
          // A := A - v*(~w) - w*(~v)
          dsyr2( kUpper, i+1, -1.0, &A(0,i+1), 1, tau, 1, pA, lda );
          A(i,i+1) = e[i];
        }
        d[i+1] = A(i+1,i+1);
        tau[i] = taui;
      }
      d[0] = A(0,0);
    }
    else
    {
      // Reduce the lower triangle of A

      for( i = 0; i < (n-1); ++i )
      {
        // Generate elementary reflector H(i) = I - tau*v*(~v)
        // to annihilate A(i+2:n-1,i)

        dlarfg( n-(i+1), A(i+1,i), &A(min(i+2,n-1),i), 1, taui );
        e[i] = A(i+1,i);

        if( 0.0 != taui )
        {
          // Apply H(i) from both sides to A(i+1:n-1,i+1:n-1)

          A(i+1,i) = 1.0;

          // Compute  x := tau*A*v  storing x in TAU(i:n-2)
          dsymv( kLower, n-(i+1), taui, &A(i+1,i+1), lda, &A(i+1,i), 1, 0.0, &tau[i], 1 );

          // Compute  w := x - (1/2 * tau*dot(x,v)) * v
          alpha = -0.5*taui*ddot( n-(i+1), &tau[i], 1, &A(i+1,i), 1 );
          daxpy( n-(i+1), alpha, &A(i+1,i), 1, &tau[i], 1 );

          // Apply the transformation as a rank-2 update:
          // A := A - v*(~w) - w*(~v)
          dsyr2( kLower, n-(i+1), -1.0, &A(i+1,i), 1, &tau[i], 1, &A(i+1,i+1), lda );
          A(i+1,i) = e[i];
        }

        d[i] = A(i,i);
        tau[i] = taui;
      }
    }
  }

  //------------------------------------------------------------------------
  // dorg2l
  //------------------------------------------------------------------------

  void dorg2l( int m, int n, int k, double *pA, int lda, double *tau, double *work )
  {
    int i, ii, j, h;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    if( m < 0 ){ throw XerblaException{ "dorg2l", 1 }; }
    if( n < 0 || n > m ){ throw XerblaException{ "dorg2l", 2 }; }
    if( k < 0 || k > n ){ throw XerblaException{ "dorg2l", 3 }; }
    if( lda < max(1,m) ){ throw XerblaException{ "dorg2l", 5 }; }

    // Quick return if possible

    if( 0 == n ){ return; }

    // Initialise columns 0:n-k to columns of the unit matrix

    for( j = 0; j < (n-k); ++j )
    {
      for( h = 0; h < m; ++h )
      { A(h,j) = 0.0; }
      A((m-n+(j+1))-1,j) = 1.0;
    }

    for( i = 0; i < k; ++i )
    {
      ii = (n-k+(i+1))-1;

      // Apply H(i) to A(0:(m-k+(i+1))-1,1:(n-k+(i+1))-1) from the left
      A(m-n+(ii+1)-1,ii) = 1.0;

      dlarf( kLeft, m-n+(ii+1), (ii-1)+1, &A(0,ii), 1, tau[i], pA, lda, work );
      dscal( m-n+(ii+1)-1, -tau[i], &A(0,ii), 1 );
      A(m-n+(ii+1)-1,ii) = 1.0 - tau[i];

      // Set A((m-k+(i+1)+1)-1:m-1,(n-k+(i+1))-1) to zero
      for( h = (m-n+(ii+1)+1)-1; h < m; ++h )
      { A(h,ii) = 0.0; }
    }
  }

  //------------------------------------------------------------------------
  // dorg2r
  //------------------------------------------------------------------------

  void dorg2r( int m, int n, int k, double *pA, int lda, double *tau, double *work )
  {
    int i, j, h;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    if( m < 0 ){ throw XerblaException{ "dorg2r", 1 }; }
    if( n < 0 || n > m ){ throw XerblaException{ "dorg2r", 2 }; }
    if( k < 0 || k > n ){ throw XerblaException{ "dorg2r", 3 }; }
    if( lda < max(1,m) ){ throw XerblaException{ "dorg2r", 5 }; }

    // Quick return if possible

    if( 0 == n ){ return; }

    // Initialise columns k+1:n to columns of the unit matrix

    for( j = (k+1)-1; j < n; ++j )
    {
      for( h = 0; h < m; ++h )
      { A(h,j) = 0.0; }
      A(j,j) = 1.0;
    }

    for( i = k-1; i >= 0; --i )
    {
      // Apply H(i) to A(i:m-1,i:n-1) from the left
      if( i < n )
      {
        A(i,i) = 1.0;
        dlarf( kLeft, m-(i+1)+1, n-(i+1), &A(i,i), 1, tau[i], &A(i,i+1), lda, work );
      }
      if( i < m )
      { dscal( m-(i+1), -tau[i], &A(i+1,i), 1 ); }

      A(i,i) = 1.0 - tau[i];

      // Set A(1:i-1,i) to zero
      for( h = 0; h <= i-1; ++h )
      { A(h,i) = 0.0; }
    }
  }

  //------------------------------------------------------------------------
  // dorgtr
  //------------------------------------------------------------------------

  // work is must be at least n-1 in size
  void dorgtr( HalfOpt uplo, int n, double *pA, int lda, double *tau, double *work )
  {
    int i, j;

    auto A = [&]( int i, int j ) noexcept -> double &
    { return pA[col_mjr(i,j,lda)]; };

    // NOTE: Removed block-based logic.

    // Test the input arguments

    if( n < 0 ){ throw XerblaException{ "dorgqr", 2 }; }

    // Quick return if possible

    if( 0 == n ){ return; }

    if( kUpper == uplo )
    {
      // Q was determined by a call to DSYTRD with UPLO = 'U'
      //
      // Shift the vectors which define the elementary reflectors one
      // column to the left, and set the last row and column of Q to
      // those of the unit matrix

      for( j = 0; j < n-1; ++j )
      {
        for( i = 0; i <= j-1; ++i )
        { A(i,j) = A(i,j+1); }
        A(n-1,j) = 0.0;
      }
      for( i = 0; i < n-1; ++i )
      { A(i,n-1) = 0.0; }
      A(n-1,n-1) = 1.0;

      // Generate Q(0:n-2,0:n-2)

      dorg2l( n-1, n-1, n-1, pA, lda, tau, work );
    }
    else if( kLower == uplo )
    {
      // Q was determined by a call to DSYTRD with UPLO = 'L'.
      //
      // Shift the vectors which define the elementary reflectors one
      // column to the right, and set the first row and column of Q to
      // those of the unit matrix

      for( j = n-1; j >= 1; --j )
      {
        A(0,j) = 0.0;
        for( i = j+1; i < n; ++i )
        { A(i,j) = A(i,j+1); }
      }
      A(0,0) = 1.0;
      for( i = 1; i < n; ++i )
      { A(i,0) = 0.0; }

      // Generate Q(1:n-1,1:n-1)

      dorg2r( n-1, n-1, n-1, &A(1,1), lda, tau, work );
    }
  }

  //------------------------------------------------------------------------
  // dsterf
  //------------------------------------------------------------------------

  // Does NOT sort the eigenvalues
  void dsterf( int n, double *d, double *e, int maxit )
  {
    // Mmmmm ...delicious goto sphagetti...
    //
    // This needed a couple of variable names changed to make it easier to read...

    int i, h, iscale, it, g, g1, gend, gendsv, hsv, nmaxit;
    double alpha, anorm, bb, c, eps, eps2, gamma, sigma, oldc, oldgam;
    double p, r, rt1, rt2, rte, s, safmax, safmin, ssfmax, ssfmin;//, rmax;

    // Test the input parameters.

    if( n < 0 ){ throw XerblaException{ "dsterf", 1 }; }

    // Quick return if possible

    if( n <= 1 ){ return; }

    eps = DBL_EPSILON;
    eps2 = sqr(eps);
    safmin = 1.0/DBL_MAX;
    safmax = DBL_MAX;
    ssfmax = sqrt( safmax )/3.0;
    ssfmin = sqrt( safmin )/eps2;
    //rmax = DBL_MAX*(1.0-DBL_EPSILON);

    // Compute the eigenvalues of the tridiagonal matrix.

    nmaxit = n*maxit;
    sigma = 0.0;
    it = 0;

    // Determine where the matrix splits and choose QL or QR iteration
    // for each block, according to whether top or bottom diagonal
    // element is smaller.

    g1 = 0;

  jmp_10__:

    if( g1 > (n-1) )
    { goto jmp_170_; }

    if( g1 > 0 )
    { e[g1-1] = 0.0; }

    for( h = g1; h < n-1; ++h )
    {
      if( fabs(e[h]) <= ( eps*sqrt(fabs(d[h]))*sqrt(fabs(d[h+1])) ) )
      {
        e[h] = 0.0;
        goto jmp_30__;
      }
    }

    h = n-1;

  jmp_30__:

    g = g1;
    hsv = g;
    gend = h;
    gendsv = gend;
    g1 = h+1;
    if( gend == g )
    { goto jmp_10__; }

    // Scale submatrix in rows and columns G to GEND

    anorm = dlanst( kMaxNorm, gend-g+1, &d[g], &e[g] );
    iscale = 0;
    if( 0.0 == anorm )
    { goto jmp_10__; }

    if( anorm > ssfmax )
    {
      iscale = 1;
      dlascl( kFullMat, 0, 0, anorm, ssfmax, gend-g+1, 1, &d[g], n );
      dlascl( kFullMat, 0, 0, anorm, ssfmax, gend-g, 1, &e[g], n );
    }
    else if( anorm < ssfmin )
    {
      iscale = 2;
      dlascl( kFullMat, 0, 0, anorm, ssfmin, gend-g+1, 1, &d[g], n );
      dlascl( kFullMat, 0, 0, anorm, ssfmin, gend-g, 1, &e[g], n );
    }

    for( i = g; i <= (gend-1); ++i )
    { e[i] *= e[i]; }

    // Choose between QL and QR iteration

    if( fabs(d[gend]) < fabs(d[g]) )
    {
      gend = hsv;
      g = gendsv;
    }

    if( gend >= g )
    {
      // QL Iteration

  jmp_50__:

      // Look for a small subdiagonal element.

      if( g != gend )
      {
        for( h = g; h <= (gend-1); ++h )
        {
          if( sqr(e[h]) <= ( eps2*fabs(d[h]*d[h+1]) + safmin ) )
          { goto jmp_70__; }
        }
      }
      h = gend;

  jmp_70__:

      if( h < gend )
      { e[h] = 0.0; }
      p = d[g];

      if( h == g )
      { goto jmp_90__; }

      // If remaining matrix is 2 by 2, use DLAE2 to compute its
      // eigenvalues.

      if( h == g+1 )
      {
        rte = sqrt(e[g]);
        dlae2( d[g], rte, d[g+1], rt1, rt2 );
        d[g] = rt1;
        d[g+1] = rt2;
        e[g] = 0.0;
        g += 2;
        if( g <= gend )
        { goto jmp_50__; }
        goto jmp_150_;
      }

      if( it == nmaxit )
      { goto jmp_150_; }
      ++it;

      // Form shift.
      rte = sqrt(e[g]);
      sigma = (d[g+1]-p)/(2.0*rte);
      r = dlapy2( sigma, 1.0 );
      sigma = p - ( rte/(sigma+copysign(r,sigma)) );
      c = 1.0;
      s = 0.0;
      gamma = d[h] - sigma;
      p = (gamma*gamma);

      // Inner loop

      for( i = h-1; i >= g; --i )
      {
        bb = e[i];
        r = p + bb;
        if( i != (h-1) )
        { e[i+1] = s*r; }
        oldc = c;
        c = p/r;
        s = bb/r;
        oldgam = gamma;
        alpha = d[i];
        gamma = c*(alpha-sigma) - s*oldgam;
        d[i+1] = oldgam + (alpha-gamma);
        if( 0.0 != c )
        { p = (gamma*gamma)/c; }else
        { p = oldc*bb; }
      }

      e[g] = s*p;
      d[g] = sigma + gamma;
      goto jmp_50__;

  jmp_90__:
      // Eigenvalue found.

      d[g] = p;

      ++g;
      if( g <= gend )
      { goto jmp_50__; }
      goto jmp_150_;
    }
    else
    {
      // QR Iteration

  jmp_100_:

      // Look for a small superdiagonal element.

      for( h = g; h >= (gend+1); --h )
      {
        if( fabs(e[h-1]) <= ( eps2*fabs(d[h]*d[h-1]) + safmin ) )
        { goto jmp_120_; }
      }
      h = gend;

  jmp_120_:

      if( h > gend )
      { e[h-1] = 0.0; }
      p = d[g];
      if( h == g )
      { goto jmp_140_; }

      // If remaining matrix is 2 by 2, use DLAE2 to compute its
      // eigenvalues.

      if( h == g-1 )
      {
        rte = sqrt(e[g-1]);
        dlae2( d[g], rte, d[g-1], rt1, rt2 );
        d[g] = rt1;
        d[g-1] = rt2;
        e[g-1] = 0.0;
        g -= 2;
        if( g >= gend )
        { goto jmp_100_; }
        goto jmp_150_;
      }

      if( it == nmaxit )
      { goto jmp_150_; }
      ++it;

      // Form shift.

      rte = sqrt(e[g-1]);
      sigma = (d[g-1]-p)/(2.0*rte);
      r = dlapy2( sigma, 1.0 );
      sigma = p - (e[g-1]/(sigma+copysign(r,sigma)));
      c = 1.0;
      s = 0.0;
      gamma = d[h] - sigma;
      p = (gamma*gamma);

      // Inner loop

      for( i = h; i <= (g-1); ++i )
      {
        bb = e[i];
        r = p + bb;
        if( i != h )
        { e[i-1] = s*r; }
        oldc = c;
        c = p/r;
        s = bb/r;
        oldgam = gamma;
        alpha = d[i+1];
        gamma = c*(alpha-sigma) - s*oldgam;
        d[i] = oldgam + (alpha-gamma);
        if( 0.0 != c )
        { p = (gamma*gamma)/c; }else
        { p = oldc*bb; }
      }

      e[g-1] = s*p;
      d[g] = sigma + gamma;
      goto jmp_100_;

  jmp_140_:
      // Eigenvalue found.

      d[g] = p;

      --g;
      if( g >= gend )
      { goto jmp_100_; }
      goto jmp_150_;
    }

    // Undo scaling if necessary

  jmp_150_:

    if( 1 == iscale )
    { dlascl( kFullMat, 0, 0, ssfmax, anorm, (gendsv+1)-(hsv+1)+1, 1, &d[hsv], n ); }
    if( 2 == iscale )
    { dlascl( kFullMat, 0, 0, ssfmin, anorm, (gendsv+1)-(hsv+1)+1, 1, &d[hsv], n ); }

    // Check for no convergence to an eigenvalue after a total
    // of N*MAXIT iterations.

    if( it < nmaxit ){ goto jmp_10__; }

    goto jmp_180_;

  jmp_170_:
  jmp_180_:

    return;
  }

  //------------------------------------------------------------------------
  // dsteqr
  //------------------------------------------------------------------------

  // Does NOT sort the eigensystem
  void dsteqr( CompEigvOpt compz, int n, double *d, double *e, double *pZ, int ldz, double *work, int maxit )
  {
    int i, iscale, jtot, l, l1, lend, lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1, nm1, nmaxit;
    double anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2, s, safmax, safmin, ssfmax, ssfmin, tst;

    // Note: All m-prefixed variables are considered in the index space, NOT the size space as otherwise implied.

    auto Z = [&]( int i, int j ) noexcept -> double &
    { return pZ[col_mjr(i,j,ldz)]; };

    // Test the input parameters.

    if( n < 0 ){ throw XerblaException{ "dsteqr", 2 }; }
    if( ldz < 1 ){ throw XerblaException{ "dsteqr", 6 }; }

    // Quick return if possible

    if( 0 == n ){ return; }

    if( ( 1 == n ) && ( kTridEigv == compz ))
    { Z(0,0) = 1.0; }

    // Determine the unit roundoff and over/underflow thresholds.

    eps = DBL_EPSILON;
    eps2 = sqr(eps);
    safmin = 1.0/DBL_MAX;
    safmax = DBL_MAX;
    ssfmax = sqrt( safmax )/3.0;
    ssfmin = sqrt( safmin )/eps2;
    //rmax = DBL_MAX*(1.0-DBL_EPSILON); <- not used

    // Compute the eigenvalues and eigenvectors of the tridiagonal
    // matrix.

    if( kTridEigv == compz )
    { dlaset( kUpAndLo, n, n, 0.0, 1.0, pZ, ldz ); }

    nmaxit = n*maxit;
    jtot = 0;

    // Determine where the matrix splits and choose QL or QR iteration
    // for each block, according to whether top or bottom diagonal
    // element is smaller.

    l1 = 0;
    nm1 = n - 1;

  jmp_10__:

    if( l1 > (n-1) )
    { goto jmp_160_; }

    if( l1 > 0 )
    { e[l1-1] = 0.0; }

    if( l1 <= (nm1-1) )
    {
      for( m = l1; m < nm1; ++m )
      {
        tst = fabs( e[m] );
        if( 0.0 == tst ){ goto jmp_30__; }
        if( tst <= ( eps*sqrt(fabs(d[m]))*sqrt(fabs(d[m+1])) ) )
        {
          e[m] = 0.0;
          goto jmp_30__;
        }
      }
    }

    m = n-1;

  jmp_30__:

    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m+1;
    if( lend == l )
    { goto jmp_10__; }

    // Scale submatrix in rows and columns L to LEND

    anorm = dlanst( kMaxNorm, lend-l+1, &d[l], &e[l] );
    iscale = 0;
    if( 0.0 == anorm )
    { goto jmp_10__; }

    if( anorm > ssfmax )
    {
      iscale = 1;
      dlascl( kFullMat, 0, 0, anorm, ssfmax, lend-l+1, 1, &d[l], n );
      dlascl( kFullMat, 0, 0, anorm, ssfmax, lend-l, 1, &e[l], n );
    }
    else if( anorm < ssfmin )
    {
      iscale = 2;
      dlascl( kFullMat, 0, 0, anorm, ssfmin, lend-l+1, 1, &d[l], n );
      dlascl( kFullMat, 0, 0, anorm, ssfmin, lend-l, 1, &e[l], n );
    }

    // Choose between QL and QR iteration

    if( fabs(d[lend]) < fabs(d[l]) )
    {
      lend = lsv;
      l = lendsv;
    }

    if( lend > l )
    {
      // QL Iteration

      // Look for a small subdiagonal element.
  jmp_40__:

      if( l != lend )
      {
        lendm1 = lend - 1;
        for( m = l; m <= lendm1; ++m )
        {
          tst = sqr(e[m]);
          if( tst <= ( eps2*fabs(d[m])*fabs(d[m+1]) ) + safmin )
          { goto jmp_60__; }
        }
      }
      m = lend;

  jmp_60__:

      if( m < lend )
      { e[m] = 0.0; }
      p = d[l];
      if( m == l )
      { goto jmp_80__; }

      // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
      // to compute its eigensystem.
      if( m == (l+1) )
      {
        if( compz > kNoEigv )
        {
          dlaev2( d[l], e[l], d[l+1], rt1, rt2, c, s );
          work[l] = c;
          work[(n-1+(l+1))-1] = s;
          dlasr( kRight, kPivotVar, kBackwd, n, 2, &work[l], &work[n-1+(l+1)-1], &Z(0,l), ldz );
        }
        else
        {
          dlae2( d[l], e[l], d[l+1], rt1, rt2 );
        }
        d[l] = rt1;
        d[l+1] = rt2;
        e[l] = 0.0;
        l += 2;
        if( l <= lend )
        { goto jmp_40__; }
        goto jmp_140_;
      }

      if( jtot == nmaxit )
      { goto jmp_140_; }
      ++jtot;

      // Form shift.
      g = ( d[l+1] - p )/( 2.0*e[l] );
      r = dlapy2( g, 1.0 );
      g = d[m] - p + ( e[l]/(g + copysign(r,g)) );
      s = 1.0;
      c = 1.0;
      p = 0.0;

      // Inner loop

      mm1 = m - 1;
      for( i = mm1; i >= l; --i )
      {
        f = s*e[i];
        b = c*e[i];
        dlartg( g, f, c, s, r );
        if( i != m-1 )
        { e[i+1] = r; }
        g = d[i+1] - p;
        r = (d[i] - g)*s + 2.0*c*b;
        p = s*r;
        d[i+1] = g + p;
        g = c*r - b;

        // If eigenvectors are desired, then save rotations.

        if( compz > kNoEigv )
        {
          work[i] = c;
          work[n-1+(i+1)-1] = -s;
        }
      }

      // If eigenvectors are desired, then apply saved rotations.

      if( compz > kNoEigv )
      {
        mm = (m-l+1)-1;
        dlasr( kRight, kPivotVar, kBackwd, n, mm+1, &work[l], &work[n-1+(l+1)-1], &Z(0,l), ldz );
      }

      d[l] -= p;
      e[l] = g;
      goto jmp_40__;

  jmp_80__:

      // Eigenvalue found.
      d[l] = p;
      ++l;
      if( l <= lend )
      { goto jmp_40__; }
      goto jmp_140_;
    }
    else
    {
      // QR Iteration

      // Look for a small superdiagonal element.
  jmp_90__:
      if( l != lend )
      {
        lendp1 = lend + 1;
        for( m = l; m >= lendp1; --m )
        {
          tst = sqr(e[m-1]);
          if( tst <= ( eps2*fabs(d[m])*fabs(d[m-1]) + safmin ) )
          { goto jmp_110_; }
        }
      }
      m = lend;

  jmp_110_:
      if( m > lend )
      { e[m-1] = 0.0; }
      p = d[l];
      if( m == l )
      { goto jmp_130_; }

      // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
      // to compute its eigensystem.

      if( m == (l-1) )
      {
        if( compz > kNoEigv )
        {
          dlaev2( d[l-1], e[l-1], d[l], rt1, rt2, c, s );
          work[m] = c;
          work[n-1+(m+1)-1] = s;
          dlasr( kRight, kPivotVar, kForwd, n, 2, &work[m], &work[n-1+(m+1)-1], &Z(0,l-1), ldz );
        }
        else
        {
          dlae2( d[l-1], e[l-1], d[l], rt1, rt2 );
        }
        d[l-1] = rt1;
        d[l] = rt2;
        e[l-1] = 0.0;
        l -= 2;
        if( l >= lend )
        { goto jmp_90__; }
        goto jmp_140_;
      }

      if( jtot == nmaxit )
      { goto jmp_140_; }
      ++jtot;

      // Form shift.

      g = ( d[l-1] - p )/( 2.0*e[l-1] );
      r = hypot( g, 1.0 );
      g = d[m] - p + ( e[l-1]/( g+copysign(r,g) ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;

      // Inner loop

      lm1 = l - 1;
      for( i = m; i <= lm1; ++i )
      {
        f = s*e[i];
        b = c*e[i];
        dlartg( g, f, c, s, r );
        if( i != m )
        { e[i-1] = r; }
        g = d[i] - p;
        r = ( d[i+1] - g )*s + 2.0*c*b;
        p = s*r;
        d[i] = g + p;
        g = c*r - b;

        // If eigenvectors are desired, then save rotations.

        if( compz > kNoEigv )
        {
          work[i] = c;
          work[n-1+(i+1)-1] = s;
        }
      }

      // If eigenvectors are desired, then apply saved rotations.

      if( compz > kNoEigv )
      {
        mm = (l-m+1)-1;
        dlasr( kRight, kPivotVar, kForwd, n, mm+1, &work[m], &work[n-1+(m+1)-1], &Z(0,m), ldz );
      }

      d[l] -= p;
      e[lm1] = g;
      goto jmp_90__;

  jmp_130_:

      // Eigenvalue found.
      d[l] = p;
      --l;
      if( l >= lend )
      { goto jmp_90__; }
      goto jmp_140_;
    }

  jmp_140_:

    // Undo scaling if necessary
    if( 1 == iscale )
    {
      dlascl( kFullMat, 0, 0, ssfmax, anorm, lendsv-lsv+1, 1, &d[lsv], n );
      dlascl( kFullMat, 0, 0, ssfmax, anorm, lendsv-lsv+1, 1, &d[lsv], n );
    }
    else if( 2 == iscale )
    {
      dlascl( kFullMat, 0, 0, ssfmin, anorm, lendsv-lsv+1, 1, &d[lsv], n );
      dlascl( kFullMat, 0, 0, ssfmin, anorm, lendsv-lsv, 1, &e[lsv], n );
    }

    // Check for no convergence to an eigenvalue after a total
    // of N*MAXIT iterations.

    if( jtot < nmaxit )
    { goto jmp_10__; }
    goto jmp_190_;

  jmp_160_:
  jmp_190_:

    return;
  }

  void TestSymmEigensolve()
  {
    double baseErrMax = 1.0e-10;

    auto CmpDbl = []( const void *pa, const void *pb ) noexcept -> int
    {
      const auto &a = *(double*)pa;
      const auto &b = *(double*)pb;
      if( a > b ){ return  1; }
      if( a < b ){ return -1; }
      return 0;
    };

    auto PrintVec = []( const char *name, int n, const double *v )
    {
      printf( "%s[] = { ", name );
      for( int i = 0; i < (n-1); ++i )
      { printf( "%g, ", v[i] ); }
      printf( "%g };\n", v[n-1] );
    };

    auto PrintSqMat = []( const char *name, int n, const double *pA )
    {
      printf( "%s[] = {\n", name );
      for( int j = 0; j < n; ++j )
      {
        printf( "{ " );
        for( int i = 0; i < (n-1); ++i )
        { printf( "%g, ", pA[col_mjr(j,i,n)] ); }
        if( j < (n-1) )
        { printf( "%g },\n", pA[col_mjr(j,n-1,n)] ); }else
        { printf( "%g }\n};\n", pA[col_mjr(j,n-1,n)] ); }
      }
    };

    auto DoTest = [&]( const char *title, int n, const double *pA0 )
    {
      std::vector<double> tmp{};
      tmp.resize( (size_t)(4*(n*n+1) + 5*(n+2) + 2*n), std::numeric_limits<double>::quiet_NaN() );

      double * pA = &tmp.front();
      double * pZ0 = pA + (n*n + 1);
      double * pZ = pZ0 + (n*n + 1);
      double * pL = pZ + (n*n + 1);
      double * d = pL + (n*n + 1);
      double * e = d + (n+2);
      double * d0 = e + (n+2);
      double * e0 = d0 + (n+2);
      double * tau = e0 + (n+2);
      double * work = tau + (n+2);
      bool failed = false;

      // two calls to gemm - multiply the error estimate by at least n...
      double errMax = baseErrMax * (n*n*n);

      printf( "RUNNING TEST: %s...\n", title );

      PrintSqMat( "A", n, pA0 );

      printf("--------\n");

      memcpy( pA, pA0, n*n*sizeof(double) );

      dsytd2( kUpper, n, pA, n, d, e, tau );

      if( ! isnan(d[n]) )
      {
        failed = true;
        printf( "FAILED: dsytd2 d[] array overrun.\n" );
      }
      if( ! isnan(e[n-1]) )
      {
        failed = true;
        printf( "FAILED: dsytd2 e[] array overrun.\n" );
      }
      if( ! isnan(pA[n*n]) )
      {
        failed = true;
        printf( "FAILED: dsytd2 A[] array overrun.\n" );
      }

      if( failed )
      { return; }

      memcpy( d0, d, n*sizeof(double) );
      memcpy( e0, e, (n-1)*sizeof(double) );

      dsterf( n, d, e, 1024 );

      if( ! isnan(d[n]) )
      {
        failed = true;
        printf( "FAILED: dsterf d[] array overrun.\n" );
      }
      if( ! isnan(e[n-1]) )
      {
        failed = true;
        printf( "FAILED: dsterf e[] array overrun.\n" );
      }

      if( failed )
      { return; }

      memcpy( pZ, pA, n*n*sizeof(double) );

      dorgtr( kUpper, n, pZ, n, tau, work );

      if( ! isnan(pZ[n*n]) )
      {
        failed = true;
        printf("FAILED: dorgtr Z[] array overrun.\n");
      }
      if( ! isnan(tau[n-1]) )
      {
        failed = true;
        printf("FAILED: dorgtr tau[] array overrun.\n");
      }

      if( failed )
      { return; }

      dsteqr( kOrigEigv, n, d0, e0, pZ, n, work , 1024 );

      if( ! isnan(d0[n]) )
      {
        failed = true;
        printf("FAILED: dsteqr d[] array overrun.\n");
      }
      if( ! isnan(e0[n-1]) )
      {
        failed = true;
        printf("FAILED: dsteqr e[] array overrun.\n");
      }

      if( failed )
      { return; }

      // Save Z for later.
      memcpy( pZ0, pZ, n*n*sizeof(double) );

      // Note: Eigenvalue ordering may be different between
      // dsterf and dsteqr.

      // Set the diagonal along L to the eigenvalues.
      memset( pL, 0, n*n*sizeof(double) );
      for( int i = 0; i < n; ++i ){ pL[col_mjr(i,i,n)] = d0[i]; }

      // Reconstruct A := Z*L*(~Z)

      // A := Z*L
      dgemm( kNoTrans, kNoTrans, n, n, n, 1.0, pZ, n, pL, n, 0.0, pA, n );

      // L = A*(~Z)
      dgemm( kNoTrans, kTrans, n, n, n, 1.0, pA, n, pZ, n, 0.0, pL, n );

      // Check that A0 ~= Z*L*(~Z)
      for( int i = 0; i < (n*n); ++i )
      {
        const auto err = fabs(pA0[i] - pL[i]);
        if( err > errMax )
        {
          failed = true;
          printf( "FAILED: mismatch at flat element %i ( | %g - %g | = %g ) in A0 and Z*L*(~Z). errMax = %g.\n",
            i, pA0[i], pL[i], err, errMax );
          break;
        }
      }

      memcpy( pL, d0, n*sizeof(double) );

      // Sort and check eigenvalues
      qsort( d, n, sizeof(double), CmpDbl );
      qsort( d0, n, sizeof(double), CmpDbl );

      for( int i = 0; i < n; ++i )
      {
        const auto err = fabs(d[i] - d0[i]);
        if( err > errMax )
        {
          failed = true;
          printf( "FAILED: eigenvalue mismatch at sorted index %i ( | %g - %g | = %g ). errMax = %g.\n",
            i, d[i], d0[i], err, errMax );
          break;
        }
      }

      if( ! failed )
      {
        printf("PASSED!\n");
        printf("Results from dsteqr...\n");
        PrintVec( "d", n, pL );
        PrintSqMat( "Z", n, pZ0 );
      }

      printf("--------\n");
    };

    double A_4x4[]
    {
      199.119, 101.671, -15.6411, -27.3375, // <- column 0
      101.671, 193.39, 84.9177, -67.7843,
      -15.6411, 84.9177, 200.791, -18.6486,
      -27.3375, -67.7843, -18.6486, 25.7863 // <- column 3
    };

    DoTest( "4x4 Random Symmetric Matrix", 4, A_4x4 );

    double A_5x5[]
    {
      1, 2, 3, 4, 5, // <- column 0
      2, 3, 4, 5, 6,
      3, 4, 5, 6, 7,
      4, 5, 6, 7, 8,
      5, 6, 7, 8, 9  // <- column 4
    };

    DoTest( "5x5 Contrived Symmetric Matrix", 5, A_5x5 );

    double A_16x16[]
    {
      5.081806042603272e6, -259190.1002349069, -1.5293764995366044e6, 602041.0408513317,
      -1.6228693880274545e6, -2.2963310390380807e6, -1.3460058149299102e6, 335589.82059287705,
      315354.6152972983, 577025.102408141, -1.2018271769795064e6, 549390.9768791221,
      -1.3130520949310476e6, -1.086981559266605e6, 1.1927771708342305e6, -769421.6059497147,
      -259190.1002349069, 6.699068391451641e6, 1.1044807246942087e6, -977492.8997798273,
      1.790195153782133e6, 75449.80160577905, 2.545644794805552e6, 1.8168743524536039e6,
      1.798478642098208e6, 75268.92931343742, 817328.6380715376, -632562.7725910144,
      -180610.83929921585, 343700.68287136895, 174795.8771230716, -2.9232664012057944e6,
      -1.5293764995366044e6, 1.1044807246942087e6, 4.753158268932328e6, -3.728268535667603e6,
      2.007736203246121e6, 710641.2875531535, -1.629787549442682e6, -882306.1244953319,
      1.2828616743240862e6, 2.240121091145197e6, -1.8738136482248972e6, -3.0833648935840963e6,
      -472132.60878733726, 1.1937187052797985e6, -721840.6018289515, -628008.2444476318,
      602041.0408513317, -977492.8997798273, -3.728268535667603e6, 6.313609454156766e6,
      -470288.7286653321, -163536.11778482073, 823811.1186597677, 1.7391625098530664e6,
      154019.73756080883, -3.536466984947088e6, 1.113990152979018e6, 1.7182253539175692e6,
      2.216875192722892e6, -1.4125180534764891e6, 237321.5372484878, 2.187163680131318e6,
      -1.6228693880274545e6, 1.790195153782133e6, 2.007736203246121e6, -470288.7286653321,
      5.933160251303321e6, 1.5369313936204428e6, -1.3378620338581514e6, -223413.67654114176,
      1.9538307374697612e6, -297074.50788362883, 1.4746604145148625e6, -2.711672912828507e6,
      -802365.208982448, -1.2485407498532431e6, 138089.2654333172, -505469.67061248166,
      -2.2963310390380807e6, 75449.80160577905, 710641.2875531535, -163536.11778482073,
      1.5369313936204428e6, 3.673097585109671e6, -258752.6119814822, -1.650741390429146e6,
      974516.759636021, 46110.56267403079, -58661.63125638288, 400912.5847754644,
      17837.631643051463, -25062.454330494467, -1.1426310909245394e6, 861805.9821513785,
      -1.3460058149299102e6, 2.545644794805552e6, -1.629787549442682e6, 823811.1186597677,
      -1.3378620338581514e6, -258752.6119814822, 6.534386882450767e6, 211683.23375603074,
      167804.25531978725, -1.38113353221188e6, 793283.0193390812, 836678.5163058355,
      763909.2698364774, 1.6355468272193277e6, 192382.48507304827, -1.0331326529332065e6,
      335589.82059287705, 1.8168743524536039e6, -882306.1244953319, 1.7391625098530664e6,
      -223413.67654114176, -1.650741390429146e6, 211683.23375603074, 4.974783667892501e6,
      -913626.4174030427, -109722.50415338005, 27826.30035456996, -409566.98049590306,
      1.278028161358026e6, -2.3577948953148685e6, 677270.7682269526, -1.3410262148926724e6,
      315354.6152972983, 1.798478642098208e6, 1.2828616743240862e6, 154019.73756080883,
      1.9538307374697612e6, 974516.759636021, 167804.25531978725, -913626.4174030427,
      5.812839086051912e6, 529002.3987502307, -737008.7820409115, -1.4190854726763438e6,
      657369.6106800145, 715629.7168073589, -710477.711607071, 738489.9951541702,
      577025.102408141, 75268.92931343742, 2.240121091145197e6, -3.536466984947088e6,
      -297074.50788362883, 46110.56267403079, -1.38113353221188e6, -109722.50415338005,
      529002.3987502307, 5.538722842640498e6, -1.5589106047466348e6, -1.2064864519751535e6,
      -309519.68957561627, 305296.489428011, 350886.2860278949, -861715.4228547735,
      -1.2018271769795064e6, 817328.6380715376, -1.8738136482248972e6, 1.113990152979018e6,
      1.4746604145148625e6, -58661.63125638288, 793283.0193390812, 27826.30035456996,
      -737008.7820409115, -1.5589106047466348e6, 4.888691404408463e6, 753107.4464271434,
      144456.00237855676, 59044.86006698321, 402814.1498676421, 971491.1458501453,
      549390.9768791221, -632562.7725910144, -3.0833648935840963e6, 1.7182253539175692e6,
      -2.711672912828507e6, 400912.5847754644, 836678.5163058355, -409566.98049590306,
      -1.4190854726763438e6, -1.2064864519751535e6, 753107.4464271434, 5.732180720806532e6,
      191072.2760868772, -467595.73363894667, -326420.7342306595, 358931.03278663,
      -1.3130520949310474e6, -180610.83929921567, -472132.60878733714, 2.216875192722892e6,
      -802365.208982448, 17837.631643051413, 763909.2698364774, 1.278028161358026e6,
      657369.6106800145, -309519.6895756162, 144456.0023785567, 191072.27608687733,
      3.213796178651998e6, -108993.03693795021, -167909.9244147583, 989511.6543502613,
      -1.0869815592666047e6, 343700.68287136895, 1.1937187052797985e6, -1.4125180534764891e6,
      -1.2485407498532434e6, -25062.454330494416, 1.635546827219328e6, -2.3577948953148685e6,
      715629.7168073588, 305296.4894280111, 59044.86006698347, -467595.7336389467,
      -108993.03693795021, 5.885608654454638e6, -311616.4284996298, 333525.51501671056,
      1.1927771708342303e6, 174795.87712307167, -721840.6018289515, 237321.53724848802,
      138089.26543331728, -1.1426310909245391e6, 192382.48507304813, 677270.7682269526,
      -710477.7116070709, 350886.2860278949, 402814.1498676421, -326420.7342306594,
      -167909.9244147583, -311616.4284996298, 1.6447614985350743e6, -393928.1477244597,
      -769421.6059497149, -2.923266401205794e6, -628008.2444476315, 2.1871636801313185e6,
      -505469.6706124815, 861805.9821513785, -1.0331326529332068e6, -1.3410262148926724e6,
      738489.9951541707, -861715.4228547731, 971491.1458501451, 358931.0327866301,
      989511.6543502613, 333525.51501671056, -393928.1477244597, 7.029275099005061e6
    };

    DoTest( "16x16 Random Symmetric Matrix", 16, A_16x16 );
  }

}; // struct Toy_Cxx_LAPACK_3_7_0
