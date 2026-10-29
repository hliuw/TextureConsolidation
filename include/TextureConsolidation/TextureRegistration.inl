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

template< unsigned int Quadrature , typename Solver >
TextureRegistration< Quadrature , Solver >::TextureRegistration( EmbeddedPhongTexturedMesh< K , K > mesh , const unsigned int res[K] , Parameters params )
	: _mesh( mesh ) , _params( params ) , _textureAdvection( mesh , res , false ) , _textureComparer( mesh )
{
	Miscellany::PerformanceMeter pMeter( '.' );
	_P = _mesh.tangentProlongation();
	_Pt = _P.transpose();
	_S = _Pt * _mesh.template stiffnessMatrix< Quadrature >() * _P;
	_flowSolver.analyzePattern( _S );
	_flowMassEntries = _mesh.eigenMatrixEntries();
	if( _params.verbose ) std::cout << pMeter( "Flow system" ) << std::endl;
}

template< unsigned int Quadrature , typename Solver >
std::vector< typename TextureRegistration< Quadrature , Solver >::FlowField > TextureRegistration< Quadrature , Solver >::operator()( const std::vector< TextureAndConfidence > &textureAndConfidences )
{
	Miscellany::PerformanceMeter pMeter( '.' );
	const unsigned int n = static_cast< unsigned int >( textureAndConfidences.size() );
	if( n<2 ) MK_ERROR_OUT( "Expected at least two textures: " , n );

	// Normalize the colors
	std::vector< TextureAndConfidence > normalized = textureAndConfidences;
	{
		TextureNormalizer< 32 , Solver , K , Dim > normalizer( _mesh , textureAndConfidences[0].texture.res(0) , textureAndConfidences[0].texture.res(1) );
		for( unsigned int i=0 ; i<n ; i++ )
		{
			normalizer.update( _params.normalizationWeight , &textureAndConfidences[i].confidence );
			normalizer( normalized[i].texture , &textureAndConfidences[i].confidence );
		}
	}
	if( _params.verbose ) std::cout << pMeter( "Normalized" ) << std::endl;

	// Initialize the flow-fields to zero and advect the normalized textures along them, eroding the confidences
	std::vector< FlowField > flows( n , FlowField( _mesh.vertexNum() ) );
	std::vector< TextureAndConfidence > flowed( n );
	for( unsigned int i=0 ; i<n ; i++ )
	{
		flowed[i] = _textureAdvection( flows[i] , normalized[i] );
		if( _params.erodeRadius ) flowed[i].confidence = _Erode( flowed[i].confidence , _params.erodeRadius , _params.erodeThreshold );
	}

	// The anchor: the blend of the normalized textures before any flow
	const bool anchored = _params.anchorWeight>0;
	if( anchored ) _anchor = StitchByWeightedConfidence( flowed );

	double blurWeight = _params.blurWeight , regularizationWeight = _params.regularizationWeight;
	double anchorBlur = _params.blurWeight , anchorWeight = _params.anchorWeight;
	for( unsigned int l=0 ; l<_params.levels ; l++ )
	{
		Miscellany::PerformanceMeter pMeter( '.' );

		// Select the target: of the pair of textures with largest overlap, the one with larger coverage
		unsigned int s=-1 , t=-1;
		{
			double maxOverlap = -std::numeric_limits< double >::infinity();
			for( unsigned int i=0 ; i<n ; i++ ) for( unsigned int j=i+1 ; j<n ; j++ )
			{
				double overlap = _textureComparer.overlap( flowed[i].confidence , flowed[j].confidence );
				if( overlap>maxOverlap ) maxOverlap = overlap , s = i , t = j;
			}
			if( _textureComparer.coverage( flowed[s].confidence )>_textureComparer.coverage( flowed[t].confidence ) ) std::swap( s , t );
		}
		double preError = _params.verbose ? _textureComparer.error( flowed ) : 0;

		std::set< unsigned int > registered = { t };
		std::string order;
		TextureAndConfidence target = flowed[t];
		while( true )
		{
			// Register the source to the target and update the flow-field of the source
			FlowField delta = anchored ? _anchoredFlow( flowed[s] , target , blurWeight , regularizationWeight , anchorBlur , anchorWeight ) : _flow( flowed[s] , target , blurWeight , regularizationWeight );
			for( size_t v=0 ; v<flows[s].size() ; v++ ) flows[s][v] += delta[v];

			// Advect the normalized source and merge it into the target
			flowed[s] = _textureAdvection( flows[s] , normalized[s] );
			registered.insert( s );
			order += " " + std::to_string( s );
			target = StitchByWeightedConfidence( flowed , registered );
			if( registered.size()==n ) break;

			// Select the next source: the one with the largest overlap with the target
			double maxOverlap = -std::numeric_limits< double >::infinity();
			for( unsigned int i=0 ; i<n ; i++ ) if( registered.find(i)==registered.end() )
			{
				double overlap = _textureComparer.overlap( target.confidence , flowed[i].confidence );
				if( overlap>maxOverlap ) maxOverlap = overlap , s = i;
			}
		}
		if( _params.verbose )
		{
			std::cout << "Level " << l << ": target " << t << ", sources" << order << ", error " << preError << " -> " << _textureComparer.error( flowed ) << std::endl;
			std::cout << pMeter( "Level " + std::to_string(l) ) << std::endl;
		}

		// Reduce the smoothing for the next level
		blurWeight *= _params.factor;
		regularizationWeight *= _params.factor;

		// The anchor weight decays by anchorDecay; the anchor smoothing follows the level's
		// smoothing for the first half of the levels and is then held, so the ghosted blend is never
		// compared at fine scales
		anchorWeight *= _params.anchorDecay;
		if( l+1<=_params.levels/2 ) anchorBlur *= _params.factor;
	}
	return flows;
}

template< unsigned int Quadrature , typename Solver >
typename TextureRegistration< Quadrature , Solver >::TextureAndConfidence TextureRegistration< Quadrature , Solver >::advect( const FlowField &flow , const TextureAndConfidence &textureAndConfidence , double minStepSize ) const
{
	return _textureAdvection( flow , textureAndConfidence , 1. , minStepSize );
}

template< unsigned int Quadrature , typename Solver >
typename TextureRegistration< Quadrature , Solver >::FlowField TextureRegistration< Quadrature , Solver >::_flow( const TextureAndConfidence &source , const TextureAndConfidence &target , double blurWeight , double regularizationWeight )
{
	// Blur the source and target, taking confidence into account, sampling them to the vertices
	TextureToVertex< Quadrature , Solver , K , Dim > sourceToVertex( _mesh ) , targetToVertex( _mesh );
	sourceToVertex.update( blurWeight , &source.confidence );
	targetToVertex.update( blurWeight , &target.confidence );
	return _flow( sourceToVertex( source.texture , &source.confidence ) , targetToVertex( target.texture , &target.confidence ) , source.confidence , target.confidence , regularizationWeight );
}

template< unsigned int Quadrature , typename Solver >
typename TextureRegistration< Quadrature , Solver >::FlowField TextureRegistration< Quadrature , Solver >::_anchoredFlow( const TextureAndConfidence &source , const TextureAndConfidence &target , double blurWeight , double regularizationWeight , double anchorBlur , double anchorWeight )
{
	// The source/target term at the level's smoothing, and the source/anchor term at the anchor smoothing
	TextureToVertex< Quadrature , Solver , K , Dim > sourceToVertex( _mesh ) , targetToVertex( _mesh ) , coarseSourceToVertex( _mesh ) , anchorToVertex( _mesh );
	sourceToVertex.update( blurWeight , &source.confidence );
	targetToVertex.update( blurWeight , &target.confidence );
	coarseSourceToVertex.update( anchorBlur , &source.confidence );
	anchorToVertex.update( anchorBlur , &_anchor.confidence );
	FitTerm anchorTerm;
	anchorTerm.sourceSignal = coarseSourceToVertex( source.texture , &source.confidence );
	anchorTerm.targetSignal = anchorToVertex( _anchor.texture , &_anchor.confidence );
	anchorTerm.sourceConfidence = &source.confidence;
	anchorTerm.targetConfidence = &_anchor.confidence;
	anchorTerm.weight = anchorWeight;
	return _flow( sourceToVertex( source.texture , &source.confidence ) , targetToVertex( target.texture , &target.confidence ) , source.confidence , target.confidence , regularizationWeight , &anchorTerm );
}

template< unsigned int Quadrature , typename Solver >
typename TextureRegistration< Quadrature , Solver >::FlowField TextureRegistration< Quadrature , Solver >::_flow
(
	const std::vector< Point< double , Channels > > &sourceSignal ,
	const std::vector< Point< double , Channels > > &targetSignal ,
	const RegularGrid< K , double > &sourceConfidence ,
	const RegularGrid< K , double > &targetConfidence ,
	double regularizationWeight ,
	const FitTerm *extra
)
{
	using T = Point< double , Channels >;

	// The flow minimizes the confidence-weighted linearized optical-flow energy:
	//		E(V) = \int w * \sum_c [ < \nabla Source[c] , V > - ( Source[c] - Target[c] ) ]^2 + beta_reg * \int || \nabla V ||^2
	// with w the product of the source and target confidences, giving the system
	//		( \int w (\nabla Source)(\nabla Source)^T + beta_reg * S ) V = \int w ( Source - Target ) \nabla Source
	// With an extra term, its fitting energy (with its own signals and confidences) is added, scaled by its weight.

	TextureScalarField< K , double , K > sourceConfidenceField( _mesh.textureCoordinates() , sourceConfidence ) , targetConfidenceField( _mesh.textureCoordinates() , targetConfidence );

	std::vector< T > signalDiff( sourceSignal.size() );
	for( size_t i=0 ; i<sourceSignal.size() ; i++ ) signalDiff[i] = sourceSignal[i] - targetSignal[i];
	ScalarField< K , T > signalDiffField( _mesh , signalDiff );
	ScalarField< K , T > sourceField( _mesh , sourceSignal );
	typename ScalarField< K , T >::Differential dSource( sourceField );

	// The same for the extra term, built only when there is one
	// (otherwise the references below alias the first term's data, and the extra term is never evaluated)
	const bool hasExtra = extra!=nullptr;
	const double extraWeight = hasExtra ? extra->weight : 0.;
	std::optional< TextureScalarField< K , double , K > > _extraSourceConfidenceField , _extraTargetConfidenceField;
	std::vector< T > extraDiff;
	std::optional< ScalarField< K , T > > _extraDiffField , _extraSourceField;
	std::optional< typename ScalarField< K , T >::Differential > _dExtraSource;
	if( hasExtra )
	{
		_extraSourceConfidenceField.emplace( _mesh.textureCoordinates() , *extra->sourceConfidence );
		_extraTargetConfidenceField.emplace( _mesh.textureCoordinates() , *extra->targetConfidence );
		extraDiff.resize( sourceSignal.size() );
		for( size_t i=0 ; i<sourceSignal.size() ; i++ ) extraDiff[i] = extra->sourceSignal[i] - extra->targetSignal[i];
		_extraDiffField.emplace( _mesh , extraDiff );
		_extraSourceField.emplace( _mesh , extra->sourceSignal );
		_dExtraSource.emplace( *_extraSourceField );
	}
	const TextureScalarField< K , double , K > &extraSourceConfidenceField = hasExtra ? *_extraSourceConfidenceField : sourceConfidenceField;
	const TextureScalarField< K , double , K > &extraTargetConfidenceField = hasExtra ? *_extraTargetConfidenceField : targetConfidenceField;
	const ScalarField< K , T > &extraDiffField = hasExtra ? *_extraDiffField : signalDiffField;
	const typename ScalarField< K , T >::Differential &dExtraSource = hasExtra ? *_dExtraSource : dSource;

	// The transformation taking intrinsic tangents to extrinsic tangents
	PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K > I2E( _mesh );

	// The weights
	auto W = [&]( size_t sIdx ){ return [s=sourceConfidenceField[sIdx],t=targetConfidenceField[sIdx]]( SimplexProcessing::Position< K > p ){ return s(p) * t(p); }; };
	auto WE = [&]( size_t sIdx ){ return [s=extraSourceConfidenceField[sIdx],t=extraTargetConfidenceField[sIdx],extraWeight]( SimplexProcessing::Position< K > p ){ return s(p) * t(p) * extraWeight; }; };

	// The outer product of the source gradients
	auto _OuterB = [&]( size_t sIdx , const typename ScalarField< K , T >::Differential &d )
		{
			return [gInv=_mesh.inverseMetricTensorField(sIdx),i2e=I2E[sIdx],ds=d[sIdx]]( SimplexProcessing::Position< K > p )
				{
					Point< Point< double , Channels > , Dim , double > g = i2e(p) * gInv(p) * ds(p);
					return OuterProduct( g , g );
				};
		};
	auto _B = [&]( size_t sIdx ){ return _OuterB( sIdx , dSource ); };
	auto _BE = [&]( size_t sIdx ){ return _OuterB( sIdx , dExtraSource ); };

	// The source gradient scaled by the difference between the source and target
	auto _GradL = [&]( size_t sIdx , const ScalarField< K , T > &diff , const typename ScalarField< K , T >::Differential &d )
		{
			return [sdf=diff[sIdx],gInv=_mesh.inverseMetricTensorField(sIdx),i2e=I2E[sIdx],ds=d[sIdx]]( SimplexProcessing::Position< K > p )
				{
					Point< Point< double , Dim > , Channels , double > g = _Transpose( i2e(p) * gInv(p) * ds(p) );
					T d = sdf(p);
					Point< double , Dim > v;
					for( unsigned int c=0 ; c<Channels ; c++ ) v += g[c] * d[c];
					return v;
				};
		};
	auto _L = [&]( size_t sIdx ){ return _GradL( sIdx , signalDiffField , dSource ); };
	auto _LE = [&]( size_t sIdx ){ return _GradL( sIdx , extraDiffField , dExtraSource ); };

	// The squared difference between the source and target
	auto _C = [&]( size_t sIdx ){ return [sdf=signalDiffField[sIdx]]( SimplexProcessing::Position< K > p ){ return sdf(p).squareNorm(); }; };
	auto _CE = [&]( size_t sIdx ){ return [sdf=extraDiffField[sIdx]]( SimplexProcessing::Position< K > p ){ return sdf(p).squareNorm(); }; };

	auto B = [&]( size_t sIdx ){ return [_b=_B(sIdx),w=W(sIdx),_be=_BE(sIdx),we=WE(sIdx),hasExtra]( SimplexProcessing::Position< K > p ){ return hasExtra ? _b(p) * w(p) + _be(p) * we(p) : _b(p) * w(p); }; };
	auto L = [&]( size_t sIdx ){ return [_l=_L(sIdx),w=W(sIdx),_le=_LE(sIdx),we=WE(sIdx),hasExtra]( SimplexProcessing::Position< K > p ){ return hasExtra ? _l(p) * w(p) + _le(p) * we(p) : _l(p) * w(p); }; };

	// The fitting matrix and dual vector
	auto SysB = [&]( size_t sIdx )
		{
			return [sys=SimplexProcessing::PhongRodriguesSystem< K , Dim >::template MassMatrix< false >( _mesh.simplexVertices(sIdx) , _mesh.simplexNormals(sIdx) , B( sIdx ) )]( SimplexProcessing::Position< K > p ){ return sys(p); };
		};
	_mesh.template setSystemMatrixEntries< Quadrature >( _flowMassEntries , SimplexProcessing::ArrayWrapper( SysB ) );
	Eigen::VectorXd flowDual = _mesh.template massVector< Quadrature >( SimplexProcessing::ArrayWrapper( L ) );

	_flowSolver.factorize( _Pt * _flowMassEntries() * _P + _S * regularizationWeight );
	if( _flowSolver.info()!=Eigen::Success ) MK_ERROR_OUT( "Failed to factorize system" );

	Eigen::VectorXd x = _P * _flowSolver.solve( _Pt * flowDual );
	if( x.size()!=_mesh.vertexNum() * Dim ) MK_ERROR_OUT( "Unexpected size: " , x.size() , " != " , _mesh.vertexNum() , " * " , Dim );
	FlowField flow( _mesh.vertexNum() );
	for( size_t i=0 ; i<flow.size() ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) flow[i][d] = x[i*Dim+d];

	// Rescale the flow by the step size minimizing the (linearized) fitting energy
	// This is the GetScale operation, see "Motion Graphs for Unstructured Textured Meshes" (Prada et al., 2016)
	{
		PhongRodriguesVectorField< K > flowField( _mesh , flow );
		auto Energy = [&]( size_t sIdx )
			{
				return [_b=_B(sIdx),_l=_L(sIdx),_c=_C(sIdx),w=W(sIdx),_be=_BE(sIdx),_le=_LE(sIdx),_ce=_CE(sIdx),we=WE(sIdx),hasExtra,ff=flowField[sIdx]]( SimplexProcessing::Position< K > p )
					{
						Polynomial::Polynomial< 1 , 2 , double > poly;
						Point< double , Dim > f = ff(p);
						poly.coefficient(2) = _b(p)(f,f);
						poly.coefficient(1) = -2 * Point< double , Dim >::Dot( f , _l(p) );
						poly.coefficient(0u) = _c(p);
						if( !hasExtra ) return poly * w(p);
						Polynomial::Polynomial< 1 , 2 , double > polyE;
						polyE.coefficient(2) = _be(p)(f,f);
						polyE.coefficient(1) = -2 * Point< double , Dim >::Dot( f , _le(p) );
						polyE.coefficient(0u) = _ce(p);
						return poly * w(p) + polyE * we(p);
					};
			};
		double scale = 1.;
		Polynomial::Polynomial1D< 2 > poly = _mesh.template integral< Quadrature , Polynomial::Polynomial1D< 2 > >( SimplexProcessing::ArrayWrapper( Energy ) );
		if( Polynomial::Roots( poly.d(0) , &scale )!=1 ) MK_ERROR_OUT( "Could not find root" );
		for( size_t i=0 ; i<flow.size() ; i++ ) flow[i] *= scale;
	}
	return flow;
}
