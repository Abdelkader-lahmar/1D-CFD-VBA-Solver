#include "patankarLimitedCorrection.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Limiter, template<class> class LimitFunc>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::patankarLimitedCorrection<Type, Limiter, LimitFunc>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "patankarLimitedCorrection::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), Zero)
        )
    );

    // tsfCorr is zero everywhere by construction. Internal faces and
    // coupled-patch faces are LEFT AT ZERO deliberately: their high-
    // order behaviour is already fully handled by weights() (inherited
    // unchanged from LimitedScheme -- the real vanLeer/QUICK/MUSCL
    // limiting logic). We only add a nonzero correction at non-coupled
    // outflow boundary faces below.

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >
    > tgradVf = gradScheme_().grad(vf, gradSchemeName_);

    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& gradVf = tgradVf();

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        Boundary& bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            // Coupled (processor/cyclic) patches: weights() already
            // handles these fully via LimitedScheme::calcLimiter's
            // coupled branch. No extra correction wanted here.
            continue;
        }

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];
        const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

        forAll(pOwner, facei)
        {
            const label own = pOwner[facei];

            if (pFaceFlux[facei] > 0)
            {
                // Outflow: extrapolate using the owner cell's
                // gradient -- same mechanism as patankarLinearUpwind.
                pSfCorr[facei] = (pCf[facei] - C[own]) & gradVf[own];
            }
            // Inflow: leave at zero. The boundary condition
            // (patankarFixedValueFvPatchField) already supplies the
            // exact prescribed value for the base convective term
            // there via valueInternalCoeffs/valueBoundaryCoeffs.
        }
    }

    return tsfCorr;
}

// ************************************************************************* //
