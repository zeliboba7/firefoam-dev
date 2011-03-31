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

#include "sequentialNew.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ODEChemistryType>
Foam::sequentialNew<ODEChemistryType>::sequentialNew
(
    const fvMesh& mesh,
    const word& ODEmodelName,
    const word& thermoType
)
:
    chemistrySolver<ODEChemistryType>(mesh, ODEmodelName, thermoType),
    coeffsDict_(this->subDict(ODEmodelName + "Coeffs")),
    cTauChem_(readScalar(coeffsDict_.lookup("cTauChem"))),
    eqRateLimiter_(coeffsDict_.lookup("equilibriumRateLimiter"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ODEChemistryType>
Foam::sequentialNew<ODEChemistryType>::~sequentialNew()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ODEChemistryType>
Foam::scalar Foam::sequentialNew<ODEChemistryType>::solve
(
    scalarField &c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{
    scalar tChemInv = SMALL;

    scalar pf, cf, pb, cb;
    label lRef, rRef;

    for (label i=0; i<this->reactions().size(); i++)
    {
        scalar om0 = this->omegaI
        (
            i, c, T, p, pf, cf, lRef, pb, cb, rRef
        );

        scalar omega = 0.0;

        if (!eqRateLimiter_)
        {
            omega = om0;
        }
        else
        {
            if (om0 < 0.0)
            {
                omega = om0/(1.0 + pb*dt);
            }
            else
            {
                omega = om0/(1.0 + pf*dt);
            }
        }

        tChemInv = max(tChemInv, mag(omega));

        this->updateConcsInReactionI(i, dt, omega, c);
    } // end for (label i...

    return cTauChem_/tChemInv;
}


// ************************************************************************* //
