#include "upstreamGrad.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * //

template<class Type>
Foam::tmp<typename Foam::fv::upstreamGrad<Type>::GradFieldType>
Foam::fv::upstreamGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<GradFieldType> tgGrad
    (
        GradFieldType::New
        (
            name,
            mesh,
            dimensioned<GradType>(vsf.dimensions()/dimLength, Zero)
        )
    );
    GradFieldType& gGrad = tgGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& C = mesh.C();

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>("phi");

    boolList done(mesh.nCells(), false);

    // Internal faces: assign the DOWNSTREAM cell's gradient using
    // (phi_downstream - phi_upstream) / distance, where "upstream" is
    // chosen from the local sign of the flux phi on that face.
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        label upCell;
        label downCell;

        if (phi[facei] >= 0)
        {
            // Flow direction: owner -> neighbour
            upCell = own;
            downCell = nei;
        }
        else
        {
            // Flow direction: neighbour -> owner
            upCell = nei;
            downCell = own;
        }

        if (!done[downCell])
        {
            const vector d = C[downCell] - C[upCell];
            const scalar dist = mag(d);

            if (dist > SMALL)
            {
                gGrad[downCell] =
                    (vsf[downCell] - vsf[upCell]) / dist * (d / dist);
                done[downCell] = true;
            }
        }
    }

    // Boundary-adjacent cells with no internal upstream neighbour yet
    // assigned (i.e. cells sitting right next to an inflow boundary):
    // use the boundary face value as the upstream point.
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];

        if (p.coupled())
        {
            continue;
        }

        const labelUList& faceCells = p.faceCells();
        const fvPatchField<Type>& psf = vsf.boundaryField()[patchi];
        const fvsPatchScalarField& phip = phi.boundaryField()[patchi];

        forAll(faceCells, i)
        {
            const label celli = faceCells[i];

            // Outward-pointing Sf convention: phi > 0 means flux is
            // LEAVING the domain at this face. phi <= 0 here means
            // the flux is ENTERING the domain, i.e. this is an inflow
            // boundary face and can supply the "upstream" value.
            if (!done[celli] && phip[i] <= 0)
            {
                const vector d = C[celli] - p.Cf()[i];
                const scalar dist = mag(d);

                if (dist > SMALL)
                {
                    gGrad[celli] =
                        (vsf[celli] - psf[i]) / dist * (d / dist);
                    done[celli] = true;
                }
            }
        }
    }

    gGrad.correctBoundaryConditions();

    return tgGrad;
}

// ************************************************************************* //
