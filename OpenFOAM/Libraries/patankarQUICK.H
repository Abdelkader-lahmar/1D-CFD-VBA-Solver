#ifndef patankarQUICK_H
#define patankarQUICK_H

#include "patankarLimitedCorrection.H"
#include "QUICK.H"
#include "NVDTVD.H"
#include "LimitFuncs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class patankarQUICK
:
    public patankarLimitedCorrection
    <
        Type,
        QUICKLimiter<NVDTVD>,
        limitFuncs::magSqr
    >
{
public:

    TypeName("patankarQUICK");

    patankarQUICK(const fvMesh& mesh, const surfaceScalarField& faceFlux)
    :
        patankarLimitedCorrection
        <
            Type, QUICKLimiter<NVDTVD>, limitFuncs::magSqr
        >
        (
            mesh, faceFlux
        )
    {}

    patankarQUICK(const fvMesh& mesh, Istream& schemeData)
    :
        patankarLimitedCorrection
        <
            Type, QUICKLimiter<NVDTVD>, limitFuncs::magSqr
        >
        (
            mesh, schemeData
        )
    {}

    patankarQUICK
    (
        const fvMesh& mesh,
        const surfaceScalarField& faceFlux,
        Istream& schemeData
    )
    :
        patankarLimitedCorrection
        <
            Type, QUICKLimiter<NVDTVD>, limitFuncs::magSqr
        >
        (
            mesh, faceFlux, schemeData
        )
    {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
