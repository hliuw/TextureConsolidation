// Part of TextureConsolidation, distributed under its BSD-3-Clause license (see LICENSE).

// Gradient-domain stitching (Sec. 6.3, Eq. 7): TextureSignalProcessing's TextureStitching,
// compiled unmodified.
//
// This file exists for one line. Misha's MK_ERROR_OUT prints a message and exits with
// status 0, so a failure inside TextureSignalProcessing (for example "Bad component" in
// Src/SimpleTriangleMesh.inl) would look like success to whoever runs the program.
// Redefining it to throw, before anything includes Misha, routes those errors to the
// try/catch in TextureStitching's main, which returns EXIT_FAILURE -- the convention
// BackProjection and TextureFlow already follow.
#define MK_ERROR_OUT( ... ) MK_THROW( __VA_ARGS__ )

#include "TextureStitching/TextureStitching.cpp"
