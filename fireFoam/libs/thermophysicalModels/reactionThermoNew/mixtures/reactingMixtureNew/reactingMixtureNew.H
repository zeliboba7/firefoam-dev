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

Class
    Foam::reactingMixtureNew

Description
    Foam::reactingMixtureNew

SourceFiles
    reactingMixtureNew.C

\*---------------------------------------------------------------------------*/

#ifndef reactingMixtureNew_H
#define reactingMixtureNew_H

#include "chemistryReader.H"
#include "multiComponentMixtureNew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class reactingMixtureNew Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class reactingMixtureNew
:
    public autoPtr<chemistryReader<ThermoType> >,
    public multiComponentMixtureNew<ThermoType>,
    public PtrList<Reaction<ThermoType> >
{
    // Private member functions

        //- Disallow default bitwise copy construct
        reactingMixtureNew(const reactingMixtureNew&);

        //- Disallow default bitwise assignment
        void operator=(const reactingMixtureNew&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        reactingMixtureNew(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~reactingMixtureNew()
    {}


    // Member functions

        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "reactingMixtureNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
