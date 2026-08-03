#ifndef patankarFixedValueFvPatchField_H
#define patankarFixedValueFvPatchField_H

#include "fixedValueFvPatchField.H"

namespace Foam
{

template<class Type>
class patankarFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
public:

    TypeName("patankarFixedValue");

    // Constructors

        patankarFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        patankarFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        patankarFixedValueFvPatchField
        (
            const patankarFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        patankarFixedValueFvPatchField
        (
            const patankarFixedValueFvPatchField<Type>&
        );

        patankarFixedValueFvPatchField
        (
            const patankarFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }

    // Member functions

        //- Return the matrix diagonal coefficients for CONVECTION (fvm::div)
        //  Reproduces phi_f = phi_P (i.e. zeroGradient-like for convection)
        virtual tmp<Field<Type>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients for CONVECTION (fvm::div)
        virtual tmp<Field<Type>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        // gradientInternalCoeffs() / gradientBoundaryCoeffs() are
        // deliberately NOT overridden: they are inherited unchanged from
        // fixedValueFvPatchField, so fvm::laplacian still uses
        // (phi_B - phi_P)/delta as the diffusive boundary treatment.
};

} // End namespace Foam

#ifdef NoRepository
    #include "patankarFixedValueFvPatchField.C"
#endif

#endif

// ************************************************************************* //
