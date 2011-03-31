/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField

Description
    Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField

SourceFiles
    uniformDensityHydrostaticTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformDensityHydrostaticTotalPressureFvPatchScalarField_H
#define uniformDensityHydrostaticTotalPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class uniformDensityHydrostaticTotalPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class uniformDensityHydrostaticTotalPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Constant density in the far-field
        scalar rhoRefValue_;

        //- Reference pressure
        scalarField pRefValue_;

        //- Reference pressure location
        vector pRefPoint_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the velocity field
        word UName_;

        //- Name of the density field
        word rhoName_;

public:

    //- Runtime type information
    TypeName("uniformDensityHydrostaticTotalPressure");


    // Constructors

        //- Construct from patch and internal field
        uniformDensityHydrostaticTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformDensityHydrostaticTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  uniformDensityHydrostaticTotalPressureFvPatchScalarField
        //  onto a new patch
        uniformDensityHydrostaticTotalPressureFvPatchScalarField
        (
            const uniformDensityHydrostaticTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformDensityHydrostaticTotalPressureFvPatchScalarField
        (
            const uniformDensityHydrostaticTotalPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformDensityHydrostaticTotalPressureFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        uniformDensityHydrostaticTotalPressureFvPatchScalarField
        (
            const uniformDensityHydrostaticTotalPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformDensityHydrostaticTotalPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the constant density in the far-field
            scalar rhoRefValue() const
            {
                return rhoRefValue_;
            }

            //- Return reference to the constant density in the far-field
            //  to allow adjustment
            scalar& rhoRefValue()
            {
                return rhoRefValue_;
            }

            //- Return the reference pressure
            scalarField pRefValue()
            {
                return pRefValue_;
            }

            //- Return reference to the reference pressure to allow adjustment
            const scalarField& pRefValue() const
            {
                return pRefValue_;
            }

            //- Return the pressure reference location
            const vector& pRefPoint() const
            {
                return pRefPoint_;
            }

            //- Return reference to the pressure reference location
            //  to allow adjustment
            vector& pRefPoint()
            {
                return pRefPoint_;
            }

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the name of the density field
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of the density field
            //  to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  using the given patch velocity field
            virtual void updateCoeffs(const vectorField& Up);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
