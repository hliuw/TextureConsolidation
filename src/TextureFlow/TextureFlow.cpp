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

// Registers a set of texture/confidence images to each other with surface optical flow
// and writes out the registered (advected) textures and confidences.

// Report errors by throwing (caught in main, which returns EXIT_FAILURE) rather than Misha's default of exiting with status 0
#define MK_ERROR_OUT( ... ) MK_THROW( __VA_ARGS__ )

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include "Include/PreProcessing.h"

#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/MultiThreading.h"
#include "Include/SimplicialMeshProcessing.h"
#include "Include/PlyIO.h"
#include "TextureConsolidation/PlanarSubdivision.h"
#include "TextureConsolidation/TexturedMesh.h"
#include "TextureConsolidation/TextureIO.h"
#include "TextureConsolidation/TextureRegistration.h"

using namespace MishaK;
using namespace MishaK::SimplicialMesh;

static const unsigned int Quadrature = 12;
using Registration = TextureRegistration< Quadrature , LDLtSolver >;
static const unsigned int K = Registration::K;
static const unsigned int Dim = Registration::Dim;

CmdLineParameter< std::string >
	Mesh( "mesh" ) ,
	InTextures( "inTextures" ) ,
	InConfidences( "inConfidences" ) ,
	OutTextures( "outTextures" ) ,
	OutConfidences( "outConfidences" );

CmdLineParameter< unsigned int >
	Levels( "levels" , 6 ) ,
	SubdivisionIterations( "sub" , 1 ) ,
	ErodeRadius( "erode" , 1 );

CmdLineParameter< double >
	BlurWeight( "sDiffusion" , 3e-4 ) ,
	RegularizationWeight( "sWeight" , 0.5 ) ,
	Factor( "factor" , 0.5 ) ,
	NormalizationWeight( "normalization" , 1e-3 ) ,
	ErodeThreshold( "erodeEps" , 0.2 ) ,
	MinStepSize( "minStep" , 1e-2 ) ,
	AnchorWeight( "anchorWeight" , 1. ) ,
	AnchorDecay( "anchorDecay" , 0.7 );

CmdLineReadable
	Verbose( "verbose" );

std::vector< CmdLineReadable * > params =
{
	&Mesh ,
	&InTextures ,
	&InConfidences ,
	&OutTextures ,
	&OutConfidences ,
	&Levels ,
	&BlurWeight ,
	&RegularizationWeight ,
	&Factor ,
	&NormalizationWeight ,
	&SubdivisionIterations ,
	&ErodeRadius ,
	&ErodeThreshold ,
	&MinStepSize ,
	&AnchorWeight ,
	&AnchorDecay ,
	&Verbose ,
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input textured mesh>\n" , Mesh.name.c_str() );
	printf( "\t --%s <input texture pattern, e.g. textures/%%d.png>\n" , InTextures.name.c_str() );
	printf( "\t --%s <input confidence pattern, e.g. textures/%%d.confidence.png>\n" , InConfidences.name.c_str() );
	printf( "\t[--%s <output texture pattern, e.g. flow/%%d.png>]\n" , OutTextures.name.c_str() );
	printf( "\t[--%s <output confidence pattern, e.g. flow/%%d.confidence.png>]\n" , OutConfidences.name.c_str() );
	printf( "\t[--%s <hierarchy levels (L)>=%d]\n" , Levels.name.c_str() , Levels.value );
	printf( "\t[--%s <smoothing of source and target (beta_blur)>=%g]\n" , BlurWeight.name.c_str() , BlurWeight.value );
	printf( "\t[--%s <smoothness of the flow (beta_reg)>=%g]\n" , RegularizationWeight.name.c_str() , RegularizationWeight.value );
	printf( "\t[--%s <per-level scaling of beta_blur and beta_reg>=%g]\n" , Factor.name.c_str() , Factor.value );
	printf( "\t[--%s <extent of the color normalization (epsilon_n)>=%g]\n" , NormalizationWeight.name.c_str() , NormalizationWeight.value );
	printf( "\t[--%s <number of 1-to-4 subdivisions of the mesh before computing the flows>=%d]\n" , SubdivisionIterations.name.c_str() , SubdivisionIterations.value );
	printf( "\t[--%s <confidence erosion radius, in texels>=%d]\n" , ErodeRadius.name.c_str() , ErodeRadius.value );
	printf( "\t[--%s <confidence erosion threshold>=%g]\n" , ErodeThreshold.name.c_str() , ErodeThreshold.value );
	printf( "\t[--%s <minimum advection step size>=%g]\n" , MinStepSize.name.c_str() , MinStepSize.value );
	printf( "\t[--%s <weight of the term anchoring the textures to the unregistered blend; higher values keep the textures closer to their original locations, possibly at the cost of weaker registration; 0 to disable>=%g]\n" , AnchorWeight.name.c_str() , AnchorWeight.value );
	printf( "\t[--%s <per-level multiplier of the anchor weight, from coarse to fine; lower values weaken the anchoring at the finer levels>=%g]\n" , AnchorDecay.name.c_str() , AnchorDecay.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

int _main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !Mesh.set || !InTextures.set || !InConfidences.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	Miscellany::PerformanceMeter pMeter( '.' );

	// Read the mesh, compute (area-weighted) vertex normals, and subdivide
	std::vector< Point< double , Dim > > vertices , normals;
	std::vector< Point< double , K > > textureCoordinates;
	std::vector< SimplexIndex< K > > simplices;
	PlyIO::ReadTexturedMesh( Mesh.value , vertices , textureCoordinates , simplices );
	normals.resize( vertices.size() );
	for( size_t t=0 ; t<simplices.size() ; t++ )
	{
		Simplex< double , Dim , K > simplex;
		for( unsigned int k=0 ; k<=K ; k++ ) simplex[k] = vertices[ simplices[t][k] ];
		Point< double , Dim > n = simplex.normal();
		for( unsigned int k=0 ; k<=K ; k++ ) normals[ simplices[t][k] ] += n;
	}
	for( size_t i=0 ; i<normals.size() ; i++ ) normals[i] /= Point< double , Dim >::Length( normals[i] );
	for( unsigned int i=0 ; i<SubdivisionIterations.value ; i++ )
	{
		PlanarSubdivider sub( simplices , vertices.size() );
		vertices = sub.vertexData( vertices );
		normals = sub.vertexData( normals );
		for( size_t j=0 ; j<normals.size() ; j++ ) normals[j] /= Point< double , Dim >::Length( normals[j] );
		textureCoordinates = sub.cornerData( textureCoordinates );
		simplices = sub();
	}
	EmbeddedPhongTexturedMesh< K , K > mesh( vertices , normals , textureCoordinates , simplices );
	if( Verbose.set ) std::cout << "Vertices/Simplices/Texture-Coordinates: " << vertices.size() << " / " << simplices.size() << " / " << textureCoordinates.size() << std::endl;

	// Read the textures and confidences
	unsigned int count = IndexedFileCount( InTextures.value );
	if( count<2 ) MK_ERROR_OUT( "Expected at least two textures: " , IndexedFileName( InTextures.value , count ) , " not found" );
	if( IndexedFileCount( InConfidences.value )!=count ) MK_ERROR_OUT( "Number of confidences does not match number of textures: " , IndexedFileCount( InConfidences.value ) , " != " , count );
	std::vector< Registration::TextureAndConfidence > textureAndConfidences( count );
	for( unsigned int i=0 ; i<count ; i++ )
	{
		textureAndConfidences[i].texture = ReadTexture( IndexedFileName( InTextures.value , i ) );
		textureAndConfidences[i].confidence = ReadConfidenceTexture( IndexedFileName( InConfidences.value , i ) );
		if( textureAndConfidences[i].texture.res(0)!=textureAndConfidences[i].confidence.res(0) || textureAndConfidences[i].texture.res(1)!=textureAndConfidences[i].confidence.res(1) ) MK_ERROR_OUT( "Texture and confidence different sized: " , i );
		if( textureAndConfidences[0].texture.res(0)!=textureAndConfidences[i].texture.res(0) || textureAndConfidences[0].texture.res(1)!=textureAndConfidences[i].texture.res(1) ) MK_ERROR_OUT( "Textures different sized: " , i );
	}
	if( Verbose.set ) std::cout << "Textures: " << count << " @ " << textureAndConfidences[0].texture.res(0) << " x " << textureAndConfidences[0].texture.res(1) << std::endl;
	if( Verbose.set ) std::cout << pMeter( "Read input" ) << std::endl;

	// Register
	Registration::Parameters registrationParams;
	registrationParams.levels = Levels.value;
	registrationParams.blurWeight = BlurWeight.value;
	registrationParams.regularizationWeight = RegularizationWeight.value;
	registrationParams.factor = Factor.value;
	registrationParams.normalizationWeight = NormalizationWeight.value;
	registrationParams.erodeRadius = ErodeRadius.value;
	registrationParams.erodeThreshold = ErodeThreshold.value;
	registrationParams.verbose = Verbose.set;
	registrationParams.anchorWeight = AnchorWeight.value;
	registrationParams.anchorDecay = AnchorDecay.value;
	Registration registration( mesh , textureAndConfidences[0].texture.res() , registrationParams );
	std::vector< Registration::FlowField > flows = registration( textureAndConfidences );
	if( Verbose.set ) std::cout << pMeter( "Registered" ) << std::endl;

	// Advect the input textures and confidences along the flows and write them out
	for( unsigned int i=0 ; i<count ; i++ )
	{
		Registration::TextureAndConfidence advected = registration.advect( flows[i] , textureAndConfidences[i] , MinStepSize.value );
		if( OutTextures.set ) WriteTexture( IndexedFileName( OutTextures.value , i ) , advected.texture );
		if( OutConfidences.set ) WriteConfidenceTexture( IndexedFileName( OutConfidences.value , i ) , advected.confidence );
	}
	if( Verbose.set ) std::cout << pMeter( "Advected" ) << std::endl;

	return EXIT_SUCCESS;
}

int main( int argc , char* argv[] )
{
	try{ return _main( argc , argv ); }
	catch( const std::exception &e )
	{
		std::cerr << "[ERROR] " << e.what() << std::endl;
		return EXIT_FAILURE;
	}
}
