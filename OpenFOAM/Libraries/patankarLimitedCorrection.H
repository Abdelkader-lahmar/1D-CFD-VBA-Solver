#ifndef patankarLimitedCorrection_H
#define patankarLimitedCorrection_H

#include "LimitedScheme.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class Limiter, template<class> class LimitFunc>
class patankarLimitedCorrection
:
    public LimitedScheme<Type, Limiter, LimitFunc>
{
    // Private Data

        word gradSchemeName_;
        tmp<fv::gradScheme<Type>> gradScheme_;


public:

    // Constructors

        patankarLimitedCorrection
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            LimitedScheme<Type, Limiter, LimitFunc>(mesh, faceFlux),
            gradSchemeName_("grad"),
            gradScheme_
            (
                new fv::gaussGrad<Type>(mesh)
            )
        {}

        //- Construct from Istream.
        //  The name of the flux field is read first (consumed by the
        //  LimitedScheme base constructor), then the name of the
        //  gradScheme to use for the boundary-outflow correction is
        //  read from the same stream afterwards.
        patankarLimitedCorrection
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            LimitedScheme<Type, Limiter, LimitFunc>(mesh, schemeData),
            gradSchemeName_(schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh,
                    mesh.gradScheme(gradSchemeName_)
                )
            )
        {}

        patankarLimitedCorrection
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            LimitedScheme<Type, Limiter, LimitFunc>
            (
                mesh, faceFlux, schemeData
            ),
            gradSchemeName_(schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh,
                    mesh.gradScheme(gradSchemeName_)
                )
            )
        {}


    // Member Functions

        virtual bool corrected() const
        {
            return true;
        }

        //- Zero everywhere except non-coupled outflow boundary faces.
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patankarLimitedCorrection.C"
#endif

#endif

// ************************************************************************* //
