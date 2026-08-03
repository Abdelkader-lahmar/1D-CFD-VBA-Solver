#include "patankarFixedValueFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"

namespace Foam
{
    typedef patankarFixedValueFvPatchField<scalar> patankarFixedValueFvPatchScalarField;

    defineNamedTemplateTypeNameAndDebug(patankarFixedValueFvPatchScalarField, 0);

    addToRunTimeSelectionTable
    (
        fvPatchScalarField,
        patankarFixedValueFvPatchScalarField,
        patch
    );

    addToRunTimeSelectionTable
    (
        fvPatchScalarField,
        patankarFixedValueFvPatchScalarField,
        patchMapper
    );

    addToRunTimeSelectionTable
    (
        fvPatchScalarField,
        patankarFixedValueFvPatchScalarField,
        dictionary
    );
}

// ************************************************************************* //
