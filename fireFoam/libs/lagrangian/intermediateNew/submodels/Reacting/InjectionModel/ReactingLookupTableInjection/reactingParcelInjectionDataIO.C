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

#include "reactingParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::reactingParcelInjectionData::reactingParcelInjectionData(Istream& is)
:
    thermoParcelInjectionData(is)
{
    is.check("reading Y's");
    is >> Y_;

    is.check("reactingParcelInjectionData(Istream& is)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const reactingParcelInjectionData& data
)
{
    os << static_cast<const thermoParcelInjectionData&>(data);

    os << data.Y_;

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, reactingParcelInjectionData& data)
{
    is >> static_cast<thermoParcelInjectionData&>(data);

    is.check("reading Y's");
    is >> data.Y_;

    is.check("operator(Istream&, reactingParcelInjectionData&)");

    return is;
}


// ************************************************************************* //
