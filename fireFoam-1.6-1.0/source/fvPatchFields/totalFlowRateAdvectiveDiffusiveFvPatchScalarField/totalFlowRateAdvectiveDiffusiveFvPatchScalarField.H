/*---------------------------------------------------------------------------*\

     This file is part of a preliminary version of the FireFOAM, LES solver for 
     fire applications, which is currently under development at FM Global and 
     OpenCFD. This file is based on OpenFOAM, the Open Source CFD Toolbox.
     
     Copyright (c) 2009 OpenCFD / FM Global        
  
 License 
       
     FireFOAM is free software; you can redistribute it and/or modify it 
     under the terms of the GNU General Public License as published by the 
     Free Software Foundation; either version 2 of the License, or (at your 
     option) any later version. 
  
     FireFOAM is distributed in the hope that it will be useful, but WITHOUT 
     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
     for more details. 
  
     You should have received a copy of the GNU General Public License 
     along with FireFOAM; if not, write to the Free Software Foundation, 
     Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA 
  
\*---------------------------------------------------------------------------*/

#ifndef totalFlowRateAdvectiveDiffusiveFvPatchScalarField_H
#define totalFlowRateAdvectiveDiffusiveFvPatchScalarField_H

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class
    totalFlowRateAdvectiveDiffusiveFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class totalFlowRateAdvectiveDiffusiveFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
    // Private data

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;


public:

    //- Runtime type information
    TypeName("totalFlowRateAdvectiveDiffusive");


    // Constructors

        //- Construct from patch and internal field
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        //  onto a new patch
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const totalFlowRateAdvectiveDiffusiveFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const totalFlowRateAdvectiveDiffusiveFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new
                totalFlowRateAdvectiveDiffusiveFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const totalFlowRateAdvectiveDiffusiveFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new
                totalFlowRateAdvectiveDiffusiveFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return reference to the name of the flux field
            word& phiName()
            {
                return phiName_;
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
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

