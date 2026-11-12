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

// Back-projects a set of renderings with known cameras onto a texture-atlased mesh,
// producing a texture image and a confidence image per rendering.

// Report errors by throwing (caught in main, which returns EXIT_FAILURE) rather than Misha's default of exiting with status 0
#define MK_ERROR_OUT( ... ) MK_THROW( __VA_ARGS__ )

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <filesystem>
#include "Include/PreProcessing.h"

#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/MultiThreading.h"
#include "Misha/Rasterizer2D.h"
#include "Misha/Texels.h"
#include "Include/PlyIO.h"
#include "TextureConsolidation/TextureIO.h"
#include "TextureConsolidation/TextureDistance.h"
#include "CameraParameters.h"

using namespace MishaK;

#include "BackProjection.inl"

const unsigned int DefaultResolution[] = { 1024 , 1024 };

CmdLineParameter< std::string >
	Mesh( "mesh" ) ,
	Cameras( "cameras" ) ,
	Images( "images" ) ,
	OutTextures( "outTextures" ) ,
	OutConfidences( "outConfidences" );

CmdLineParameterArray< unsigned int , 2 >
	Resolution( "res" , DefaultResolution );

CmdLineParameter< unsigned int >
	DepthRadius( "depthRadius" , 20 );

CmdLineParameter< double >
	Cutoff( "cutoff" , 0.05 );

CmdLineReadable
	Verbose( "verbose" );

std::vector< CmdLineReadable * > params =
{
	&Mesh ,
	&Cameras ,
	&Images ,
	&OutTextures ,
	&OutConfidences ,
	&Resolution ,
	&DepthRadius ,
	&Cutoff ,
	&Verbose ,
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input textured mesh>\n" , Mesh.name.c_str() );
	printf( "\t --%s <input cameras (.json)>\n" , Cameras.name.c_str() );
	printf( "\t[--%s <directory the camera image paths are relative to>=<directory of the cameras file>]\n" , Images.name.c_str() );
	printf( "\t[--%s <output texture pattern, e.g. textures/%%d.png>]\n" , OutTextures.name.c_str() );
	printf( "\t[--%s <output confidence pattern, e.g. textures/%%d.confidence.png>]\n" , OutConfidences.name.c_str() );
	printf( "\t[--%s <texture resolution>=%d %d]\n" , Resolution.name.c_str() , Resolution.values[0] , Resolution.values[1] );
	printf( "\t[--%s <depth-discontinuity falloff radius, in pixels>=%d]\n" , DepthRadius.name.c_str() , DepthRadius.value );
	printf( "\t[--%s <confidences below this are set to zero>=%g]\n" , Cutoff.name.c_str() , Cutoff.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

int _main( int argc , char* argv[] )
{
	using T = Point< double , 3 >;

	CmdLineParse( argc-1 , argv+1 , params );
	if( !Mesh.set || !Cameras.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	Miscellany::PerformanceMeter pMeter( '.' );

	BackProjection::Mesh mesh;
	PlyIO::ReadTexturedMesh( Mesh.value , mesh.vertices , mesh.tCoordinates , mesh.triangles );
	if( Verbose.set ) std::cout << "Vertices/Triangles: " << mesh.vertices.size() << " / " << mesh.triangles.size() << std::endl;

	std::vector< BackProjection::Rendering< T > > renderings;
	{
		std::filesystem::path imageDirectory = Images.set ? std::filesystem::path( Images.value ) : std::filesystem::path( Cameras.value ).parent_path();
		std::vector< CameraParameters > cameraParameters = CameraParameters::Read( Cameras.value );
		renderings.resize( cameraParameters.size() );
		for( unsigned int i=0 ; i<cameraParameters.size() ; i++ )
		{
			renderings[i].cParams = cameraParameters[i];
			std::string imageFile = ( imageDirectory / cameraParameters[i].path ).string();
			if( Verbose.set ) std::cout << "Reading image: " << imageFile << std::endl;
			// Images are stored top-to-bottom, the grid bottom-to-top
			renderings[i].image = BackProjection::FlipVertical( ReadTexture( imageFile ) );
		}
	}
	if( Verbose.set ) std::cout << pMeter( "Read input" ) << std::endl;

	BackProjection::Parameters params;
	params.res[0] = Resolution.values[0];
	params.res[1] = Resolution.values[1];
	params.depthRadius = DepthRadius.value;
	params.cutoff = Cutoff.value;
	params.verbose = Verbose.set;
	std::vector< BackProjection::TextureMap< T > > textureMaps = BackProjection::Process( mesh , renderings , params );

	for( unsigned int i=0 ; i<textureMaps.size() ; i++ )
	{
		if( OutTextures.set ) WriteTexture( IndexedFileName( OutTextures.value , i ) , textureMaps[i].texture );
		if( OutConfidences.set ) WriteConfidenceTexture( IndexedFileName( OutConfidences.value , i ) , textureMaps[i].confidence );
	}
	if( Verbose.set ) std::cout << pMeter( "Back-projected" ) << std::endl;

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
