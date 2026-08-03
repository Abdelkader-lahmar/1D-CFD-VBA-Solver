#include "patankarVanLeer.H"
#include "patankarQUICK.H"
#include "patankarMUSCL.H"
#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makelimitedSurfaceInterpolationTypeScheme(patankarVanLeer, scalar)
    makelimitedSurfaceInterpolationTypeScheme(patankarQUICK, scalar)
    makelimitedSurfaceInterpolationTypeScheme(patankarMUSCL, scalar)
}

// ************************************************************************* //
