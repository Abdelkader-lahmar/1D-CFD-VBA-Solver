#include "patankarLinearUpwind.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::patankarLinearUpwind<Type>::correction
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
                "patankarLinearUpwind::correction(" + vf.name() + ')',
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

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

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

    // Internal faces: unchanged from stock linearUpwind. Extrapolate
    // from whichever cell is upstream of the face.
    forAll(faceFlux, facei)
    {
        label celli = (faceFlux[facei] > 0) ? owner[facei] : neighbour[facei];
        sfCorr[facei] = (Cf[facei] - C[celli]) & gradVf[celli];
    }

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        Boundary& bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];
        const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

        if (pSfCorr.coupled())
        {
            // Unchanged from stock linearUpwind: coupled (processor/
            // cyclic) patches need the neighbour-side gradient for
            // inflow faces, since the "upstream" cell then lives on
            // the other side of the coupled boundary.
            const Field<typename outerProduct<vector, Type>::type>
                pGradVfNei
                (
                    gradVf.boundaryField()[patchi].patchNeighbourField()
                );

            vectorField pd(mesh.boundary()[patchi].delta());

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] = (pCf[facei] - C[own]) & gradVf[own];
                }
                else
                {
                    pSfCorr[facei] =
                        (pCf[facei] - pd[facei] - C[own]) & pGradVfNei[facei];
                }
            }
        }
        else
        {
            // THE FIX: non-coupled (physical) boundary patches.
            //
            // Stock linearUpwind leaves pSfCorr at Zero here
            // unconditionally, meaning outflow boundary faces are
            // always first-order accurate regardless of gradient
            // scheme.
            //
            // Here, for OUTFLOW faces (pFaceFlux > 0), we apply the
            // same extrapolation used for internal/coupled-outflow
            // faces: it only needs the LOCAL owner cell's gradient,
            // which is available at any physical patch, coupled or
            // not.
            //
            // For INFLOW faces (pFaceFlux <= 0), we deliberately do
            // NOT extrapolate: at a fixedValue-type inlet the
            // prescribed value is already exact, and there is no
            // upstream interior cell to extrapolate from anyway.
            // pSfCorr stays at its initialised value of Zero for
            // these faces.
            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] = (pCf[facei] - C[own]) & gradVf[own];
                }
            }
        }
    }

    return tsfCorr;
}

// ************************************************************************* //
