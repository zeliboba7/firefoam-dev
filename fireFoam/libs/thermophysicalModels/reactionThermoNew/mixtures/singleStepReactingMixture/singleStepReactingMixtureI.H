/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "singleStepReactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ThermoType>
inline const Foam::dimensionedScalar
Foam::singleStepReactingMixture<ThermoType>::stoicRatio() const
{
    return stoicRatio_;
}


template<class ThermoType>
inline const Foam::dimensionedScalar
Foam::singleStepReactingMixture<ThermoType>::s() const
{
    return s_;
}


template<class ThermoType>
inline const Foam::dimensionedScalar
Foam::singleStepReactingMixture<ThermoType>::qFuel() const
{
    return qFuel_;
}


template<class ThermoType>
inline const Foam::List<Foam::scalar>&
Foam::singleStepReactingMixture<ThermoType>::specieStoichCoeffs() const
{
    return specieStoichCoeffs_;
}


template<class ThermoType>
inline Foam::tmp<Foam::volScalarField>
Foam::singleStepReactingMixture<ThermoType>::fres
(
    const label index
) const
{
    return fres_[index];
}


template<class ThermoType>
inline Foam::label
Foam::singleStepReactingMixture<ThermoType>::inertIndex() const
{
    return inertIndex_;
}


template<class ThermoType>
inline Foam::label
Foam::singleStepReactingMixture<ThermoType>::fuelIndex() const
{
    return fuelIndex_;
}


template<class ThermoType>
inline const Foam::List<int>&
Foam::singleStepReactingMixture<ThermoType>::specieProd() const
{
    return specieProd_;
}


// ************************************************************************* //
