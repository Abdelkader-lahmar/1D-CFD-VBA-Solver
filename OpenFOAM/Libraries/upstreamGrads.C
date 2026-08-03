#include "upstreamGrad.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    // ---- scalar ----
    typedef upstreamGrad<scalar> upstreamGradScalar;

    defineTemplateTypeNameAndDebugWithName
    (
        upstreamGradScalar,
        "upstream",
        0
    );

    typedef gradScheme<scalar> gradSchemeScalar;

    gradSchemeScalar::addIstreamConstructorToTable<upstreamGradScalar>
        addupstreamGradScalarIstreamConstructorToTable_;

    // ---- vector ----
    typedef upstreamGrad<vector> upstreamGradVector;

    defineTemplateTypeNameAndDebugWithName
    (
        upstreamGradVector,
        "upstream",
        0
    );

    typedef gradScheme<vector> gradSchemeVector;

    gradSchemeVector::addIstreamConstructorToTable<upstreamGradVector>
        addupstreamGradVectorIstreamConstructorToTable_;
}
}

// ************************************************************************* //
