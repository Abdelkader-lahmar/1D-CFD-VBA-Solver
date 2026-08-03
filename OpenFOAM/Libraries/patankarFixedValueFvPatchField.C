#include "patankarFixedValueFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::patankarFixedValueFvPatchField<Type>::patankarFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}

template<class Type>
Foam::patankarFixedValueFvPatchField<Type>::patankarFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict)
{}

template<class Type>
Foam::patankarFixedValueFvPatchField<Type>::patankarFixedValueFvPatchField
(
    const patankarFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}

template<class Type>
Foam::patankarFixedValueFvPatchField<Type>::patankarFixedValueFvPatchField
(
    const patankarFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf)
{}

template<class Type>
Foam::patankarFixedValueFvPatchField<Type>::patankarFixedValueFvPatchField
(
    const patankarFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patankarFixedValueFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // Same as zeroGradientFvPatchField: internal coeff = 1
    return tmp<Field<Type>>::New(this->size(), pTraits<Type>::one);
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patankarFixedValueFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // Same as zeroGradientFvPatchField: boundary source coeff = 0
    return tmp<Field<Type>>::New(this->size(), Zero);
}
