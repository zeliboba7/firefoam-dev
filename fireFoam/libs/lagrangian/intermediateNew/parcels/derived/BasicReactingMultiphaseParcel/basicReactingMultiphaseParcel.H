/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

Class
    Foam::basicReactingMultiphaseParcel

Description


SourceFiles
    basicReactingMultiphaseParcel.C

\*---------------------------------------------------------------------------*/

#ifndef basicReactingMultiphaseParcel_H
#define basicReactingMultiphaseParcel_H

#include "ReactingMultiphaseParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class basicReactingMultiphaseParcel Declaration
\*---------------------------------------------------------------------------*/

class basicReactingMultiphaseParcel
:
    public ReactingMultiphaseParcel<basicReactingMultiphaseParcel>
{

public:

    //- Run-time type information
    TypeName("basicReactingMultiphaseParcel");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        basicReactingMultiphaseParcel
        (
             ReactingMultiphaseCloud<basicReactingMultiphaseParcel>& owner,
             const vector& position,
             const label cellI
        );

        //- Construct from components
        basicReactingMultiphaseParcel
        (
             ReactingMultiphaseCloud<basicReactingMultiphaseParcel>& owner,
             const vector& position,
             const label cellI,
             const label typeId,
             const scalar nParticle0,
             const scalar d0,
             const vector& U0,
             const scalarField& YGas0,
             const scalarField& YLiquid0,
             const scalarField& YSolid0,
             const scalarField& Y0,
             const ReactingMultiphaseParcel<basicReactingMultiphaseParcel>::
                constantProperties& constProps
        );

        //- Construct from Istream
        basicReactingMultiphaseParcel
        (
            const Cloud<basicReactingMultiphaseParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        basicReactingMultiphaseParcel(const basicReactingMultiphaseParcel& p);

        //- Construct and return a clone
        autoPtr<basicReactingMultiphaseParcel> clone() const
        {
            return
                autoPtr<basicReactingMultiphaseParcel>
                (
                    new basicReactingMultiphaseParcel(*this)
                );
        }


    //- Destructor
    virtual ~basicReactingMultiphaseParcel();
};


template<>
inline bool contiguous<basicReactingMultiphaseParcel>()
{
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
