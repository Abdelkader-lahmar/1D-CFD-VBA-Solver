#ifndef upstreamGrad_H
#define upstreamGrad_H

#include "gradScheme.H"

namespace Foam
{
namespace fv
{

template<class Type>
class upstreamGrad
:
    public fv::gradScheme<Type>
{
public:

    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;


    //- Runtime type information
    TypeName("upstream");


    // Constructors

        //- Construct from mesh
        upstreamGrad(const fvMesh& mesh, Istream&)
        :
            fv::gradScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the gradient of the given field, calculated using a
        //  one-sided (backward, upstream) difference based on the local
        //  sign of the flux field "phi".
        virtual tmp<GradFieldType> calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "upstreamGrad.C"
#endif

#endif

// ************************************************************************* //
