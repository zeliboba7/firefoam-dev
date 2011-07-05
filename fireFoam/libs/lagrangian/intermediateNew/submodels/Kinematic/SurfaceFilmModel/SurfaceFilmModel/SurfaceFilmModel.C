/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "SurfaceFilmModel.H"
#include "mathematicalConstants.H"
#include "mapDistribute.H"
#include "surfaceFilmModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    g_(dimensionedVector("zero", dimAcceleration, vector::zero)),
    coeffDict_(dictionary::null),
    ejectedParcelType_(0),
    injectorCellsPatch_(0),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    UFilmPatch_(0),
    rhoFilmPatch_(0),
    deltaFilmPatch_(0),
    nParcelsTransferred_(0),
    nParcelsInjected_(0)
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const dictionary& dict,
    CloudType& owner,
    const dimensionedVector& g,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    g_(g),
    coeffDict_(dict.subDict(type + "Coeffs")),
    ejectedParcelType_(coeffDict_.lookupOrDefault("ejectedParcelType", -1)),
    injectorCellsPatch_(0),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    UFilmPatch_(0),
    rhoFilmPatch_(0),
    deltaFilmPatch_(owner.mesh().boundary().size()),
    nParcelsTransferred_(0),
    nParcelsInjected_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::~SurfaceFilmModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
template<class TrackData>
void Foam::SurfaceFilmModel<CloudType>::inject(TrackData& td)
{
    typedef regionModels::surfaceFilmModels::surfaceFilmModel filmModelType;

    bool modelPresent =
        this->owner().db().objectRegistry::foundObject<filmModelType>
        (
            "surfaceFilmProperties"
        );

    if (!modelPresent)
    {
        return;
    }

    // Retrieve the film model from the owner database
    const filmModelType& filmModel =
        this->owner().db().objectRegistry::lookupObject<filmModelType>
        (
            "surfaceFilmProperties"
        );

    if (!filmModel.active())
    {
        return;
    }

    const labelList& filmPatches = filmModel.intCoupledPatchIDs();
    const labelList& primaryPatches = filmModel.primaryPatchIDs();

    const fvMesh& mesh = owner_.mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(filmPatches, i)
    {
        const label filmPatchI = filmPatches[i];
        const label primaryPatchI = primaryPatches[i];

        const directMappedPatchBase& mapPatch =
            filmModel.mappedPatches()[filmPatchI];
        const mapDistribute& distMap = mapPatch.map();

        injectorCellsPatch_ = pbm[primaryPatchI].faceCells();

        cacheFilmFields(filmPatchI, primaryPatchI, distMap, filmModel);

        const vectorField& Cf = mesh.C().boundaryField()[primaryPatchI];
        const vectorField& Sf = mesh.Sf().boundaryField()[primaryPatchI];
        const scalarField& magSf = mesh.magSf().boundaryField()[primaryPatchI];

        forAll(injectorCellsPatch_, j)
        {
            if (diameterParcelPatch_[j] > 0)
            {
//                const point& pos = this->owner().mesh().C()[cellI];
                const scalar offset =
                    max
                    (
                        diameterParcelPatch_[j],
                        deltaFilmPatch_[primaryPatchI][j]
                    );
                const point pos = Cf[j] - 1.1*offset*Sf[j]/magSf[j];
                const label cellI = injectorCellsPatch_[j];

                // Create a new parcel
                typename CloudType::parcelType* pPtr =
                    new typename CloudType::parcelType(td.cloud(), pos, cellI);
                setParcelProperties(*pPtr, j);

                // Check new parcel properties
//                td.cloud().checkParcelProperties(*pPtr, 0.0, true);
                td.cloud().checkParcelProperties(*pPtr, 0.0, false);

                // Add the new parcel to the cloud
                td.cloud().addParticle(pPtr);

                nParcelsInjected_++;
            }
        }
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::cacheFilmFields
(
    const label filmPatchI,
    const label primaryPatchI,
    const mapDistribute& distMap,
    const regionModels::surfaceFilmModels::surfaceFilmModel& filmModel
)
{
    massParcelPatch_ = filmModel.cloudMassTrans().boundaryField()[filmPatchI];
    distMap.distribute(massParcelPatch_);

    diameterParcelPatch_ =
        filmModel.cloudDiameterTrans().boundaryField()[filmPatchI];
    distMap.distribute(diameterParcelPatch_);

    UFilmPatch_ = filmModel.U().boundaryField()[filmPatchI];
    distMap.distribute(UFilmPatch_);

    rhoFilmPatch_ = filmModel.rho().boundaryField()[filmPatchI];
    distMap.distribute(rhoFilmPatch_);

    deltaFilmPatch_[primaryPatchI] =
        filmModel.delta().boundaryField()[filmPatchI];
    distMap.distribute(deltaFilmPatch_[primaryPatchI]);
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFaceI
) const
{
    const scalar pi = mathematicalConstant::pi;

    // Set parcel properties
    scalar vol = pi/6.0*pow3(diameterParcelPatch_[filmFaceI]);
    p.d() = diameterParcelPatch_[filmFaceI];
    p.U() = UFilmPatch_[filmFaceI];
    p.rho() = rhoFilmPatch_[filmFaceI];

    p.nParticle() = massParcelPatch_[filmFaceI]/p.rho()/vol;

    if (ejectedParcelType_ >= 0)
    {
        p.typeId() = ejectedParcelType_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewSurfaceFilmModel.C"

// ************************************************************************* //
