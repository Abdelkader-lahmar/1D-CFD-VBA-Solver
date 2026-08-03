#ifndef patankarVanLeer_H
#define patankarVanLeer_H

#include "patankarLimitedCorrection.H"
#include "vanLeer.H"
#include "NVDTVD.H"
#include "LimitFuncs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class patankarVanLeer
:
    public patankarLimitedCorrection
    <
        Type,
        vanLeerLimiter<NVDTVD>,
        limitFuncs::magSqr
    >
{
public:

    TypeName("patankarVanLeer");

    patankarVanLeer(const fvMesh& mesh, const surfaceScalarField& faceFlux)
    :
        patankarLimitedCorrection
        <
            Type, vanLeerLimiter<NVDTVD>, limitFuncs::magSqr
        >
        (
            mesh, faceFlux
        )
    {}

    patankarVanLeer(const fvMesh& mesh, Istream& schemeData)
    :
        patankarLimitedCorrection
        <
            Type, vanLeerLimiter<NVDTVD>, limitFuncs::magSqr
        >
        (
            mesh, schemeData
        )
    {}

    patankarVanLeer
    (
        const fvMesh& mesh,
        const surfaceScalarField& faceFlux,
        Istream& schemeData
    )
    :
        patankarLimitedCorrection
        <
            Type, vanLeerLimiter<NVDTVD>, limitFuncs::magSqr
        >
        (
            mesh, faceFlux, schemeData
        )
    {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
