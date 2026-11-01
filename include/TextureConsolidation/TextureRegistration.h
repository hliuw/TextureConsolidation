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

#ifndef TEXTURE_REGISTRATION_INCLUDED
#define TEXTURE_REGISTRATION_INCLUDED

#include <vector>
#include <optional>
#include <set>
#include <iostream>
#include <Eigen/Sparse>
#include "Misha/Miscellany.h"
#include "Misha/Geometry.h"
#include "Misha/RegularGrid.h"
#include "Misha/Polynomial.h"
#include "Include/SimplicialMeshProcessing.h"
#include "TexturedMesh.h"
#include "TextureDistance.h"
#include "TextureToVertex.h"
#include "TextureNormalization.h"
#include "TextureAdvection.h"
#include "TextureAndConfidence.h"
#include "TextureComparison.h"

namespace MishaK
{
	namespace SimplicialMesh
	{
		// Registers a set of texture/confidence images to each other with surface optical flow (Alg. 1 of the paper):
		// the textures are normalized, then registered hierarchically, with each level successively registering
		// every texture to a composite target that is built up from the textures registered so far.
		// The result is a per-vertex (extrinsic, tangent) flow-field for each texture.
		template< unsigned int Quadrature , typename Solver >
		struct TextureRegistration
		{
			static const unsigned int K = 2 , Dim = K+1 , Channels = 3;
			using TextureAndConfidence = SimplicialMesh::TextureAndConfidence< K , Point< double , Channels > >;
			using FlowField = std::vector< Point< double , Dim > >;

			struct Parameters
			{
				unsigned int levels = 6;				// L: number of hierarchy levels
				double blurWeight = 3e-4;				// beta_blur: extent of the smoothing of source and target
				double regularizationWeight = 0.5;		// beta_reg: smoothness of the correcting flow
				double factor = 0.5;					// per-level scaling of beta_blur and beta_reg
				double normalizationWeight = 1e-3;		// epsilon_n: extent of the local mean/scale estimate used for normalization
				unsigned int erodeRadius = 1;			// confidences are eroded by this many texels before registration...
				double erodeThreshold = 0.2;			// ... from the boundary of the region where they exceed this value
				bool verbose = false;

				// Anchoring to the (unregistered) blend, so the consolidated texture does not drift: every
				// source-to-target solve also fits the source to the anchor, the confidence-weighted blend of the
				// normalized textures before any flow. The blend is ghosted, so only its low frequencies are
				// trusted: the anchor term uses the level's smoothing for the first half of the levels and then
				// holds it, and its weight starts at anchorWeight and is scaled by anchorDecay at every level.
				double anchorWeight = 1.;				// lambda_a: weight of the anchor term at level 0 (0 disables anchoring)
				double anchorDecay = 0.7;				// d: per-level multiplier of the anchor weight
			};

			TextureRegistration( EmbeddedPhongTexturedMesh< K , K > mesh , const unsigned int res[K] , Parameters params );

			// Returns the flow-field registering each of the textures
			std::vector< FlowField > operator()( const std::vector< TextureAndConfidence > &textureAndConfidences );

			// Advects a texture and its confidence along a flow-field
			TextureAndConfidence advect( const FlowField &flow , const TextureAndConfidence &textureAndConfidence , double minStepSize ) const;

		protected:
			EmbeddedPhongTexturedMesh< K , K > _mesh;
			Parameters _params;
			TextureAdvection< K > _textureAdvection;
			TextureComparer< Quadrature , K , Dim > _textureComparer;
			Eigen::SparseMatrix< double > _P , _Pt;	// Prolongation from the per-vertex tangent space to the embedding space, and its transpose
			Eigen::SparseMatrix< double > _S;		// The (restricted) covariant stiffness, regularizing the flow
			Solver _flowSolver;
			typename EmbeddedPhongMesh< K >::EigenMatrixEntries _flowMassEntries;

			// Smooths the source and target textures to the vertices and returns the flow registering the source to the target
			FlowField _flow( const TextureAndConfidence &source , const TextureAndConfidence &target , double blurWeight , double regularizationWeight );

			// A second fitting term: the source, smoothed to the vertices, against a fixed reference
			struct FitTerm
			{
				std::vector< Point< double , Channels > > sourceSignal , targetSignal;
				const RegularGrid< K , double > *sourceConfidence = nullptr , *targetConfidence = nullptr;
				double weight = 0;
			};

			// The anchor: the blend of the unflowed normalized textures (set in operator() when anchoring)
			TextureAndConfidence _anchor;

			// As _flow above, with an additional fit of the source to the anchor, smoothed with anchorBlur and weighted by anchorWeight
			FlowField _anchoredFlow( const TextureAndConfidence &source , const TextureAndConfidence &target , double blurWeight , double regularizationWeight , double anchorBlur , double anchorWeight );

			// Solves for the flow, weighted by the product of the source and target confidences, and rescales it by a line search.
			// If extra is given, its (weighted) fitting energy is added to that of the source and target.
			FlowField _flow
			(
				const std::vector< Point< double , Channels > > &sourceSignal ,
				const std::vector< Point< double , Channels > > &targetSignal ,
				const RegularGrid< K , double > &sourceConfidence ,
				const RegularGrid< K , double > &targetConfidence ,
				double regularizationWeight ,
				const FitTerm *extra = nullptr
			);

			// Zeros the confidence within the given radius (in texels) of the boundary of the region where it exceeds epsilon
			static RegularGrid< K , double > _Erode( const RegularGrid< K , double > &confidence , unsigned int radius , double epsilon )
			{
				RegularGrid< K , double > eroded = confidence;
				unsigned int radius2 = radius * radius;
				SquareDistanceToTextureBoundary< K > d2b( eroded.res() );
				RegularGrid< K , unsigned int > d2 = d2b( [&]( typename RegularGrid< K >::Index I ){ return eroded( I )>epsilon; } );
				ThreadPool::ParallelFor( 0 , eroded.size() , [&]( size_t i ){ if( d2[i]<radius2 ) eroded[i] = 0; } );
				return eroded;
			}

			static Point< Point< double , Dim > , Channels , double > _Transpose( const Point< Point< double , Channels > , Dim , double > g )
			{
				Point< Point< double , Dim > , Channels , double > t;
				for( unsigned int d=0 ; d<Dim ; d++ ) for( unsigned int c=0 ; c<Channels ; c++ ) t[c][d] = g[d][c];
				return t;
			}
		};

#include "TextureRegistration.inl"
	}
}
#endif // TEXTURE_REGISTRATION_INCLUDED
