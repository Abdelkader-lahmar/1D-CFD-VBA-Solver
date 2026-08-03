#include "upstreamGrad.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    typedef upstreamGrad<scalar> upstreamGradScalar;

    defineTemplateTypeNameAndDebugWithName
    (
        upstreamGradScalar,
        "upstream",
        0
    );

    // gradScheme's run-time selection table key is "Istream" (see
    // gradScheme.H: declareRunTimeSelectionTable(tmp, gradScheme,
    // Istream, (const fvMesh&, Istream&), (mesh, schemeData)) ), so the
    // adder member is addIstreamConstructorToTable, NOT addMeshConstructorToTable.
    // The base class must be a bare typedef (no <...>) so the member-template
    // call parses correctly.
    typedef gradScheme<scalar> gradSchemeScalar;

    gradSchemeScalar::addIstreamConstructorToTable<upstreamGradScalar>
        addupstreamGradScalarIstreamConstructorToTable_;
}
}

// ************************************************************************* //
