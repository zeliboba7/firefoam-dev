/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::interpolationLookUpTable

Description
    A list of lists. Interpolates based on the first dimension.
    The values must be positive and monotonically increasing in each dimension

    The handling of out-of-bounds values depends on the current setting
    of @a outOfBounds.

    If @a REPEAT is chosen for the out-of-bounds handling, the final time
    value is treated as being equivalent to time=0 for the following periods.

Note
    - Accessing an empty list results in an error.
    - Accessing a list with a single element always returns the same value.

SourceFiles
    interpolationLookUpTable.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationLookUpTable_H
#define interpolationLookUpTable_H

#include "List.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interpolationLookUpTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationLookUpTable
:

	  public List<List<Type> >
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum boundsHandling
        {
            ERROR,          /*!< Exit with a FatalError */
            WARN,           /*!< Issue warning and clamp value (default) */
            CLAMP,          /*!< Clamp value to the start/end value */
        };

	

private:

    // Private data

        //- Enumeration for handling out-of-bound values
        boundsHandling boundsHandling_;
        //- File name
        fileName fileName_;
        
	//- table dimensions
        List<label>  dim_;
        
	//- min on each dimension
        List<scalar>  min_;

	//- deltas on each dimension
		List<scalar>  delta_;
	
	//- maximum on each dimension
        List<scalar>  max_;

        //- entries of the Look Up Table
        List<dictionary> entries_;
		
		//- output from the Look Up Table
		List<dictionary> output_;
		
		
		HashTable<label> fieldIndices_;
		
        // Private Member Functions
        
		//- Read the table of data from file
        void readTable();

		//- Dimension Table from dictionaries input and output
		void dimensionTable();

        //- Index table finding
		label index(const List<scalar>&) const;
		
		//- Index field finding
		inline label findField(const word& field) const
		{
			return fieldIndices_[field];
		}
public:

    // Constructors

        //- Construct null
        interpolationLookUpTable();

        //- Construct given the name of the file containing the table of data
        interpolationLookUpTable(const fileName& fn);

         //- Construct from dictionary
        interpolationLookUpTable(const dictionary& dict);

        //- Construct copy
        interpolationLookUpTable(const interpolationLookUpTable& interpTable);

    // Member Functions

		// - Insert output and input fields
		void insertFields(List<dictionary>&, List<dictionary>&);

        //- Return the out-of-bounds handling as a word
        word boundsHandlingToWord(const boundsHandling& bound) const;

        //- Return the out-of-bounds handling as an enumeration
        boundsHandling wordToBoundsHandling(const word& bound) const;

        //- Set the out-of-bounds handling from enum, return previous setting
        boundsHandling outOfBounds(const boundsHandling& bound);

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void check() const;
	
        //- Write Look Up Table to filename.
        void write(Ostream& os) const;
         
		
      // Acces

		 inline const List<dictionary>& output() const
        {
            return output_;
        }

		inline const List<dictionary>& entries() const
        {
            return entries_;
        }

        inline const List<scalar>& min() const
        {
            return min_;
        }
        
        inline const List<label>& dim() const
        {
            return dim_;
        }
        
		inline const List<scalar>& delta() const
		{
	    	return delta_;
		}
	
        inline const List<scalar>& max() const
        {
            return max_;
        }
		
     // Member Operators
        
        //- Return an element of constant List<scalar, Type>
        const List<Type> & operator[](const label) const;
		
		//- Return an element of List<scalar, Type>
		List<Type> & operator[](const label);
		
        //- Return an output field  from an input vector.
		scalar operator()(scalarList& vector, const word&) const;
		
		//- Return the i output given a single input scalar (T)
		inline scalar operator()(const scalar , const label) const;
		
		//- Return the i output given two input scalars (cft, ftVar) in FireMixture
		inline scalar operator()(const scalar ,const scalar, const label) const;
		
		//- Return reference to the position count on the field "field"
		inline scalar& operator()(label count, const word& ifield)
		{
			label i = findField(ifield);
			return List<List<scalar> >::operator[](count).operator[](i);
		}

		//- operators for the output fields.
        interpolationLookUpTable<Type>& operator*=
		(
			const scalar
		);
		
		interpolationLookUpTable<Type> & operator=
		(
			const interpolationLookUpTable<Type>&
		);

		interpolationLookUpTable<Type>& operator+=
		(
			const interpolationLookUpTable&
		);

		interpolationLookUpTable<Type>& operator-=
		(
			const interpolationLookUpTable&
		);
		
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "interpolationLookUpTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
