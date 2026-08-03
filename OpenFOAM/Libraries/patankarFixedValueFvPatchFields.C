#include "patankarFixedValueFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"

namespace Foam
{
    // ---- scalar ----
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

    // ---- vector ----
    typedef patankarFixedValueFvPatchField<vector> patankarFixedValueFvPatchVectorField;

    defineNamedTemplateTypeNameAndDebug(patankarFixedValueFvPatchVectorField, 0);

    addToRunTimeSelectionTable
    (
        fvPatchVectorField,
        patankarFixedValueFvPatchVectorField,
        patch
    );

    addToRunTimeSelectionTable
    (
        fvPatchVectorField,
        patankarFixedValueFvPatchVectorField,
        patchMapper
    );

    addToRunTimeSelectionTable
    (
        fvPatchVectorField,
        patankarFixedValueFvPatchVectorField,
        dictionary
    );
}

// ************************************************************************* //
