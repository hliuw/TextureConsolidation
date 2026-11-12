/*
Copyright (c) 2026, Michael Kazhdan and Hongyi Liu
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

// Back-projection of renderings onto a texture-atlased mesh.
// For each rendering, produces a texture image (colors sampled from the rendering) and a confidence image
// (normal alignment times depth-discontinuity falloff, zero for occluded texels, quantized below a cutoff).

namespace BackProjection
{
	struct Mesh
	{
		std::vector< Point< double , 3 > > vertices;
		std::vector< Point< double , 2 > > tCoordinates;
		std::vector< SimplexIndex< 2 > > triangles;
	};

	template< typename T >
	struct Rendering
	{
		CameraParameters cParams;
		RegularGrid< 2 , T > image;
	};

	template< typename T >
	struct TextureMap
	{
		RegularGrid< 2 , double > confidence;
		RegularGrid< 2 , T > texture;
	};

	struct Parameters
	{
		unsigned int res[2];
		unsigned int depthRadius;	// r_dc: confidence falls off within this many pixels of a depth discontinuity
		double depthPrecision;		// relative tolerance of the occlusion test
		double cutoff;				// confidences below this are set to zero
		bool verbose;
		Parameters( void ) : depthRadius(20) , depthPrecision(1e-3) , cutoff(0.05) , verbose(false) { res[0] = res[1] = 1024; }
	};

	template< typename T >
	RegularGrid< 2 , T > FlipVertical( const RegularGrid< 2 , T > & g )
	{
		RegularGrid< 2 , T >  _g( g.res() );
		for( unsigned int j=0 ; j<g.res(1) ; j++ ) for( unsigned int i=0 ; i<g.res(0) ; i++ ) _g( i , j ) = g( i , g.res(1)-1-j );
		return _g;
	}

	// Rasterizes the mesh into the camera, returning the depth at each pixel (infinity where nothing is visible)
	template< typename T >
	RegularGrid< 2 , double > GetDepth( const Rendering< T > & rendering , const Mesh & mesh )
	{
		// Since we're sampling this as an unshifted RegularGrid, values are at the corners.
		static const bool NodeAtCellCenter = false;

		RegularGrid< 2 , double > depth( rendering.image.res() );
		for( size_t i=0 ; i<depth.size() ; i++ ) depth[i] = std::numeric_limits< double >::infinity();

		typename RegularGrid< 2 >::Range range;
		for( unsigned int d=0 ; d<2 ; d++ ) range.second[d] = rendering.image.res(d) - 1;

		XForm< double , 3 > kInv = rendering.cParams.k.inverse();

		for( unsigned int t=0 ; t<mesh.triangles.size() ; t++ )
		{
			Simplex< double , 2 , 2 > tTri;
			Simplex< double , 3 , 2 > cTri;
			for( unsigned int k=0 ; k<3 ; k++ )
			{
				tTri[k] = rendering.cParams( mesh.vertices[ mesh.triangles[t][k] ] );
				cTri[k] = rendering.cParams.w2c( mesh.vertices[ mesh.triangles[t][k] ] );
			}

			auto Kernel = [&]( typename RegularGrid< 2 >::Index I )
				{
					Point< double , 2 > n = kInv( Texels< NodeAtCellCenter >::NodePosition( I ) );
					Ray< double , 3 > ray;
					ray.direction = Point< double , 3 >( n[0] , n[1] , 1 );
					Point< double , 3 > bc = cTri.barycentricCoordinates( ray ).second;
					Point< double , 3 > p_c = cTri[0] * bc[0] + cTri[1] * bc[1] + cTri[2] * bc[2];
					double d = p_c[2];
					if( d>0 && d<depth(I) ) depth(I) = d;
				};
			Rasterizer2D::RasterizeNodes< NodeAtCellCenter >( tTri , Kernel , range );
		}
		return depth;
	}

	// A bilinearly-interpolated depth map, returning infinity if any of the contributing pixels is empty
	struct DepthMap
	{
		DepthMap( RegularGrid< 2 , double > && depthMap ) : _depthMap( std::move( depthMap ) ){}

		double operator()( Point< double , 2 > p ) const
		{
			double v = 0;

			int ix = static_cast< int >( floor( p[0] ) ) , iy = static_cast< int >( floor( p[1] ) );
			double dx[] = { 1. - ( p[0] - ix ) , p[0] - ix } , dy[] = { 1.0 - ( p[1] - iy ) , p[1] - iy };
			for( int ii=0 ; ii<2 ; ii++ ) for( int jj=0 ; jj<2 ; jj++ )
			{
				if( (ix+ii)>=0 && (ix+ii)<static_cast< int >( _depthMap.res(0) ) && (iy+jj)>=0 && (iy+jj)<static_cast< int >( _depthMap.res(1) ) )
				{
					if( _depthMap(ix+ii,iy+jj)<std::numeric_limits< double >::infinity() ) v += _depthMap(ix+ii,iy+jj) * dx[ii] * dy[jj];
					else return std::numeric_limits< double >::infinity();
				}
			}
			return v;
		};
	protected:
		RegularGrid< 2 , double > _depthMap;
	};

	// The depth confidence: min( 1 , D_disc / r_dc ), with D_disc the distance to the nearest depth discontinuity.
	// Discontinuities are pixels where the depth Laplacian, scaled by its root-mean-square, exceeds one.
	inline RegularGrid< 2 , double > GetDepthConfidence( const RegularGrid< 2 , double > &depth , unsigned int radius )
	{
		RegularGrid< 2 , double > laplacian( depth.res() );
		{
			typename RegularGrid< 2 >::Range range;
			for( unsigned int k=0 ; k<2 ; k++ ) range.second[k] = depth.res(k);

			auto Kernel = [&]( unsigned int , typename RegularGrid< 2 >::Index I )
				{
					typename RegularGrid< 2 >::Range nbrs = typename RegularGrid< 2 >::Range( I ).dilate( 1 );

					unsigned int count = 0;
					double value = 0;
					RegularGrid< 2 >::Range::Intersect( nbrs , range ).process( [&]( typename RegularGrid< 2 >::Index J ){ if( I!=J ){ count++ , value += depth(J); } } );
					laplacian( I ) = fabs( depth(I) - value/count );
				};
			range.processParallel( Kernel );
		}

		// Reduce in pixel order: atomic floating-point additions are race-free but
		// schedule-dependent, and the resulting threshold controls the confidence mask.
		double var = 0;
		size_t count = 0;
		for( size_t i=0 ; i<laplacian.size() ; i++ )
			if( laplacian[i]<std::numeric_limits< double >::infinity() ) count++ , var += laplacian[i] * laplacian[i];
		double dev = sqrt( var / count );
		ThreadPool::ParallelFor( 0 , laplacian.size() , [&]( size_t i ){ laplacian[i] /= dev; } );

		RegularGrid< 2 , double > depthConfidence( depth.res() );
		unsigned int radius2 = radius * radius;
		SquareDistanceToTextureBoundary< 2 > d2b( depthConfidence.res() );
		RegularGrid< 2 , unsigned int > d2 = d2b( [&]( typename RegularGrid< 2 >::Index I ){ return laplacian( I )<1; } );
		ThreadPool::ParallelFor
		(
			0 , depthConfidence.size() ,
			[&]( size_t i )
			{
				if( d2[i]<radius2 ) depthConfidence[i] = sqrt( static_cast< double >(d2[i])/radius2 );
				else                depthConfidence[i] = 1.;
			}
		);
		return depthConfidence;
	}

	template< typename T >
	std::vector< TextureMap< T > > Process( const Mesh &mesh , const std::vector< Rendering< T > > &renderings , Parameters params )
	{
		static const bool NodeAtCellCenter = true;
		using TexelInfo = typename Texels< NodeAtCellCenter >::template TexelInfo< 2 >;

		std::vector< TextureMap< T > > textureMaps( renderings.size() );

		Miscellany::PerformanceMeter pMeter( '.' );

		// For each texel, the closest point on the mesh (the triangle and barycentric coordinates)
		auto GetTextureTriangle = [&]( size_t t )
			{
				Simplex< double , 2 , 2 > tri;
				for( unsigned int k=0 ; k<3 ; k++ ) tri[k] = mesh.tCoordinates[t*3+k];
				return tri;
			};
		RegularGrid< 2 , TexelInfo > tInfo = Texels< NodeAtCellCenter >::template GetSupportedTexelInfo< 3 , false >
			(
				mesh.triangles.size() ,
				[&]( size_t v ){ return mesh.vertices[v]; } ,
				[&]( size_t t ){ return mesh.triangles[t]; } ,
				GetTextureTriangle ,
				params.res ,
				0 ,
				params.verbose
			);
		if( params.verbose ) std::cout << pMeter( "Texel info" ) << std::endl;

		for( unsigned int c=0 ; c<renderings.size() ; c++ )
		{
			RegularGrid< 2 , T > &texture = textureMaps[c].texture;
			RegularGrid< 2 , double > &confidence = textureMaps[c].confidence;
			texture.resize( params.res );
			confidence.resize( params.res );
			for( size_t i=0 ; i<confidence.size() ; i++ ) confidence[i] = 0;

			RegularGrid< 2 , double > depth = GetDepth( renderings[c] , mesh );
			RegularGrid< 2 , double > depthConfidence = GetDepthConfidence( depth , params.depthRadius );
			DepthMap dMap( std::move( depth ) );

			Point< double , 3 > cameraPosition = renderings[c].cParams.position();

			ThreadPool::ParallelFor
			(
				0 , tInfo.size() ,
				[&]( size_t i )
				{
					if( tInfo[i].sIdx!=static_cast< size_t >(-1) )
					{
						SimplexIndex< 2 > si = mesh.triangles[ tInfo[i].sIdx ];
						Simplex< double , 3  , 2 > s( mesh.vertices[ si[0] ] , mesh.vertices[ si[1] ] , mesh.vertices[ si[2] ] );
						Point< double , 3 > n = s.normal();
						Point< double , 3 > p_w = s( tInfo[i].bc );
						Point< double , 3 > dir = p_w - cameraPosition;
						Point< double , 2 > q = renderings[c].cParams( p_w );
						Point< double , 3 > p_c = renderings[c].cParams.w2c( p_w );
						double d = dMap( q );
						n /= Point< double , 3 >::Length( n );
						dir /= Point< double , 3 >::Length( dir );

						// Visible texels get the normal alignment times the depth confidence; occluded ones get zero
						if( d<std::numeric_limits< double >::infinity() && p_c[2]<d*(1.+params.depthPrecision) ) confidence[i] = fabs( Point< double , 3 >::Dot( n , dir ) ) * depthConfidence( q );
						else confidence[i] = 0.;
						if( confidence[i]<params.cutoff ) confidence[i] = 0.;

						texture[i] = confidence[i] ? renderings[c].image(q) : T();
					}
				}
			);
			if( params.verbose ) std::cout << pMeter( "Rendering " + std::to_string(c) ) << std::endl;
		}
		return textureMaps;
	}
}
