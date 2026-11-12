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

template< unsigned int K >
TextureAdvection< K >::TextureAdvection
(
	EmbeddedPhongTexturedMesh< K , K > mesh ,
	const unsigned int res[K] ,
	bool verbose
)
	: _mesh(mesh) , _rMesh( GetPointer( const_cast< std::vector< SimplexIndex< K > > & >( _mesh.simplices() ) ) , _mesh.simplexNum() ) , verbose(verbose)
{
	Miscellany::PerformanceMeter pMeter( '.' );
	_rMesh.template setMetricFromEmbedding< K+1 >( [&]( unsigned int v ){ return _mesh.vertex(v); } );
	_xForms = _rMesh.getCoordinateXForms();
	if( verbose ) std::cout << pMeter( "Riemannian mesh" ) << std::endl;

	auto TextureSimplexFunctor = [&]( size_t sIdx )
		{
			Simplex< double , K , K > simplex;
			for( unsigned int k=0 ; k<=K ; k++ ) simplex[k] = mesh.textureCoordinate( sIdx*(K+1) + k );
			return simplex;
		};
	_texelInfo = Texels< true >::template GetSupportedTexelInfo< K+1 , false >( _mesh.simplexNum() , [&]( size_t v ){ return _mesh.vertex(v); } , [&]( size_t s ){ return _mesh.simplex(s); } , TextureSimplexFunctor , res , 0 , false );
	if( verbose ) std::cout << pMeter( "Texel info" ) << std::endl;
}

template< unsigned int K >
TextureAdvection< K >::~TextureAdvection
(
	void
)
{
	DeletePointer( _xForms );
}

template< unsigned int K >
template< typename T , HasTextureAdvectionGetter< K , T > Getter , HasTextureAdvectionSetter< T > Setter >
void
TextureAdvection< K >::_process
(
	const std::vector< Point< double , K+1 > > &flowField ,
	Getter && getF ,
	Setter && setF ,
	Parameters params
)
const
{
	_TangentVectorField tvf( _mesh , flowField );

	ThreadPool::ParallelFor
	(
		0 , _texelInfo.resolution() ,
		[&]( unsigned int t , size_t i )
		{
			if( _texelInfo[i].sIdx!=-1 )
			{
				FEM::SamplePoint< double > sample;
				sample.tIdx = static_cast< int >( _texelInfo[i].sIdx );
				for( unsigned int k=0 ; k<K ; k++ ) sample.p[k] = _texelInfo[i].bc[k+1];
				_rMesh.flow( _xForms , tvf , -params.time , sample , params.minStepSize );
				setF( t , i , getF( sample.tIdx , sample.p ) );
			}
		}
	);
}

template< unsigned int K >
template< typename T >
RegularGrid< K , T >
TextureAdvection< K >::operator()
(
	const std::vector< Point< double , K+1 > > & flowField ,
	const RegularGrid< K , T > &texture ,
	Parameters params
)
const
{
	if( !params.time ) return texture;

	RegularGrid< K , T > out( texture.res() );
	// _process only writes texels supported by the atlas. RegularGrid leaves scalar
	// storage uninitialized, so confidence outside the atlas must be zero before
	// interpolation or erosion can read it (Point-valued colors already default to zero).
	for( size_t i=0 ; i<out.size() ; i++ ) out[i] = T{};
	TextureScalarField< K , T , K > textureField( _mesh.textureCoordinates() , texture );
	_process< T >( flowField , [&]( size_t s , Point< double , K > p ){ return textureField[s](p); } , [&]( unsigned int , size_t i , T t ){ out[i] = t; } , params );
	return out;
}

template< unsigned int K >
template< typename T >
TextureAndConfidence < K, T >
TextureAdvection< K >::operator()
(
	const std::vector< Point< double, K + 1 > >& flow,
	const TextureAndConfidence < K, T >& input,
	double time,
	double minStepSize
	)
	const
{
	TextureAndConfidence < K, T > out = input;
	Parameters params( minStepSize , time );
	out.texture = operator()(flow, input.texture, params);
	out.confidence = operator()(flow, input.confidence, params);
	return out;
}

template< unsigned int K >
template< typename ... Ts >
void
TextureAdvection< K >::process
(
	const std::vector< Point< double , K+1 > > &flowField ,
	const RegularGrid< K , Ts >& ... in ,
	RegularGrid< K , Ts >& ... out ,
	Parameters params
)
const
{
	_Copy< Ts... >( in... , out... );
	if( !params.time ) return;
	else
	{
		unsigned int res[K];
		_SetRes< true , Ts... >( in... , out... , res );

		auto Get = [&]( size_t sIdx , Point< double , K > bc )
			{
				Simplex< double , K , K > s;
				for( unsigned int k=0 ; k<=K ; k++ )
				{
					s[k] = _mesh.textureCoordinate( sIdx*(K+1) + k );
					for( unsigned int _k=0 ; _k<K ; _k++ ) s[k][_k] *= res[_k];
				}
				Point< double , K > p = s[0];
				for( unsigned int k=0 ; k<K ; k++ ) p += ( s[k+1]-s[0] ) * bc[k];

				std::tuple< Ts... > values;
				_ToTuple< 0 >( p , values , in... );

				return values;
			};
		auto Set = [&]( unsigned int , size_t i , std::tuple< Ts... > values ){ _FromTuple< 0 >( i , values , out... ); };

		_process< std::tuple< Ts... > >( flowField , Get , Set , params );
	}
}

template< unsigned int K >
template< unsigned int I , typename ... _T , typename T , typename ...Ts >
void TextureAdvection< K >::_ToTuple
(
	Point< double , K > p ,
	std::tuple< _T... > &values ,
	const RegularGrid< K , T > &t ,
	const RegularGrid< K , Ts >& ... ts
)
{
	static_assert( I<sizeof...(_T) , "[ERROR] Bad indexing" );
#ifdef CLAMPED_EVALUATION
#ifdef FORCE_BILINEAR
	std::get< I >( values ) = t( p );
#else // !FORCE_BILINEAR
	std::get< I >( values ) = t( p , true );
#endif // FORCE_BILINEAR
#else // !CLAMPED_EVALUATION
	std::get< I >( values ) = t( p , true ).value();
#endif // CLAMPED_EVALUATION
	if constexpr( sizeof...(Ts) ) _ToTuple< I+1 >( p , values , ts... );
}

template< unsigned int K >
template< unsigned int I , typename ... _T , typename T , typename ...Ts >
void TextureAdvection< K >::_FromTuple
(
	size_t i ,
	const std::tuple< _T... > &values ,
	RegularGrid< K , T > &t ,
	RegularGrid< K , Ts >& ... ts
)
{
	static_assert( I<sizeof...(_T) , "[ERROR] Bad indexing" );
	t[i] = std::get< I >( values );
	if constexpr( sizeof...(Ts) ) _FromTuple< I+1 >( i , values , ts... );
}

template< unsigned int K >
template< typename T , typename ... Ts >
void TextureAdvection< K >::_Copy
(
	const RegularGrid< K , T > &in ,
	const RegularGrid< K , Ts >& ... ins ,
	RegularGrid< K , T > &out ,
	RegularGrid< K , Ts >& ... outs
)
{
	out = in;
	if constexpr( sizeof...(Ts) ) _Copy< Ts... >( ins... , outs... );
}

template< unsigned int K >
template< bool First , typename T , typename ... Ts >
void TextureAdvection< K >::_SetRes
(
	const RegularGrid< K , T > &in ,
	const RegularGrid< K , Ts >& ... ins ,
	RegularGrid< K , T > &out ,
	RegularGrid< K , Ts >& ... outs , 
	unsigned int res[K]
)
{
	if constexpr( First ) for( unsigned int k=0 ; k<K ; k++ ) res[k] = in.res(k);
	else for( unsigned int k=0 ; k<K ; k++ ) if( res[k]!=in.res(k) ) MK_ERROR_OUT( "Resolution mismatch" );
	out.resize( res );
	if constexpr( sizeof...(Ts) ) _SetRes< false , Ts... >( ins... , outs... , res );
}
