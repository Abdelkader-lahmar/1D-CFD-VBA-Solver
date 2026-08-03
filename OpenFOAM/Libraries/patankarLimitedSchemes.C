#include "patankarVanLeer.H"
#include "patankarQUICK.H"
#include "patankarMUSCL.H"
#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makelimitedSurfaceInterpolationTypeScheme(patankarVanLeer, scalar)
    makelimitedSurfaceInterpolationTypeScheme(patankarVanLeer, vector)

    makelimitedSurfaceInterpolationTypeScheme(patankarQUICK, scalar)
    makelimitedSurfaceInterpolationTypeScheme(patankarQUICK, vector)

    makelimitedSurfaceInterpolationTypeScheme(patankarMUSCL, scalar)
    makelimitedSurfaceInterpolationTypeScheme(patankarMUSCL, vector)
}

// ************************************************************************* //
