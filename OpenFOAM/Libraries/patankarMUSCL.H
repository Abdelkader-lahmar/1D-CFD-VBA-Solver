#ifndef patankarMUSCL_H
#define patankarMUSCL_H

#include "patankarLimitedCorrection.H"
#include "MUSCL.H"
#include "NVDTVD.H"
#include "LimitFuncs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class patankarMUSCL
:
    public patankarLimitedCorrection
    <
        Type,
        MUSCLLimiter<NVDTVD>,
        limitFuncs::magSqr
    >
{
public:

    TypeName("patankarMUSCL");

    patankarMUSCL(const fvMesh& mesh, const surfaceScalarField& faceFlux)
    :
        patankarLimitedCorrection
        <
            Type, MUSCLLimiter<NVDTVD>, limitFuncs::magSqr
        >
        (
            mesh, faceFlux
        )
    {}

    patankarMUSCL(const fvMesh& mesh, Istream& schemeData)
    :
        patankarLimitedCorrection
        <
            Type, MUSCLLimiter<NVDTVD>, limitFuncs::magSqr
        >
        (
            mesh, schemeData
        )
    {}

    patankarMUSCL
    (
        const fvMesh& mesh,
        const surfaceScalarField& faceFlux,
        Istream& schemeData
    )
    :
        patankarLimitedCorrection
        <
            Type, MUSCLLimiter<NVDTVD>, limitFuncs::magSqr
        >
        (
            mesh, faceFlux, schemeData
        )
    {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
