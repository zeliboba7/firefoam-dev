/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "RanzMarshall.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::RanzMarshall<CloudType>::RanzMarshall
(
    const dictionary& dict,
    CloudType& cloud
)
:
    HeatTransferModel<CloudType>(dict, cloud, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::RanzMarshall<CloudType>::~RanzMarshall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
bool Foam::RanzMarshall<CloudType>::active() const
{
    return true;
}


template <class CloudType>
Foam::scalar Foam::RanzMarshall<CloudType>::Nu
(
    const scalar Re,
    const scalar Pr
) const
{
    return 2.0 + 0.6*sqrt(Re)*cbrt(Pr);
}


// ************************************************************************* //
