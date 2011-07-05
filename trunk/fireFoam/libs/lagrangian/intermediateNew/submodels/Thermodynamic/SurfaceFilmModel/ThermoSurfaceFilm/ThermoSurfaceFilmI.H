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

#include "ThermoSurfaceFilm.H"
#include "DimensionedFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::word&
Foam::ThermoSurfaceFilm<CloudType>::filmRegionName() const
{
    return filmRegionName_;
}


template<class CloudType>
inline const Foam::wordList&
Foam::ThermoSurfaceFilm<CloudType>::patchNames() const
{
    return patchNames_;
}


template<class CloudType>
inline const Foam::polyMesh&
Foam::ThermoSurfaceFilm<CloudType>::filmRegion() const
{
    return filmRegion_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::ThermoSurfaceFilm<CloudType>::hf() const
{
    return hf_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::ThermoSurfaceFilm<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField&
Foam::ThermoSurfaceFilm<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::ThermoSurfaceFilm<CloudType>::p() const
{
    return p_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::ThermoSurfaceFilm<CloudType>::h() const
{
    return h_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::ThermoSurfaceFilm<CloudType>::mu() const
{
    return mu_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoSurfaceFilm<CloudType>::rhoSp()
{
    return rhoSp_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::ThermoSurfaceFilm<CloudType>::USp()
{
    return USp_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoSurfaceFilm<CloudType>::hSp()
{
    return hSp_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoSurfaceFilm<CloudType>::pSp()
{
    return pSp_;
}


// ************************************************************************* //
