#ifndef patankarLinearUpwind_H
#define patankarLinearUpwind_H

#include "upwind.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class patankarLinearUpwind
:
    public upwind<Type>
{
    // Private Data

        word gradSchemeName_;
        tmp<fv::gradScheme<Type>> gradScheme_;


    // Private Member Functions

        patankarLinearUpwind(const patankarLinearUpwind&);
        void operator=(const patankarLinearUpwind&);


public:

    TypeName("patankarLinearUpwind");


    // Constructors

        patankarLinearUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            upwind<Type>(mesh, faceFlux),
            gradSchemeName_("grad"),
            gradScheme_
            (
                new fv::gaussGrad<Type>(mesh)
            )
        {}

        patankarLinearUpwind
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, schemeData),
            gradSchemeName_(schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh,
                    mesh.gradScheme(gradSchemeName_)
                )
            )
        {
            if (!schemeData.eof())
            {
                IOWarningInFunction(schemeData)
                    << "unexpected additional entries in stream." << nl
                    << "    Only the name of the gradient scheme in the"
                       " 'gradSchemes' dictionary should be specified."
                    << endl;
            }
        }

        patankarLinearUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, faceFlux, schemeData),
            gradSchemeName_(schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh,
                    mesh.gradScheme(gradSchemeName_)
                )
            )
        {
            if (!schemeData.eof())
            {
                IOWarningInFunction(schemeData)
                    << "unexpected additional entries in stream." << nl
                    << "    Only the name of the gradient scheme in the"
                       " 'gradSchemes' dictionary should be specified."
                    << endl;
            }
        }


    // Member Functions

        virtual bool corrected() const
        {
            return true;
        }

        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};

} // End namespace Foam

#ifdef NoRepository
    #include "patankarLinearUpwind.C"
#endif

#endif

// ************************************************************************* //
