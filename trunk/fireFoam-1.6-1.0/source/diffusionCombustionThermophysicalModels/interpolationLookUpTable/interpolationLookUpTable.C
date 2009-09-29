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

\*---------------------------------------------------------------------------*/

#include "IFstream.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
template<class Type>
Foam::label Foam::interpolationLookUpTable<Type>::index
(
  const List<scalar>& indices
) const
{
    label totalindex = 0;
    
    for(int i = 0;i < dim_.size()-1;i++)
    {
        label dim = 1;
        for(int j = i + 1 ;j <  dim_.size() ; j++)
		{
			dim *=dim_[j];
		}
        totalindex += Foam::max(label((indices[i]-min_[i])/delta_[i])-1, 0)*dim;
    }
		
    label iLastdim = dim_.size() -1;
    
    totalindex  += Foam::max
	(
		label((indices[iLastdim]-min_[iLastdim])/delta_[iLastdim])-1,
		0
	);
	
    return totalindex;
}

template<class Type>
void Foam::interpolationLookUpTable<Type>::dimensionTable()
{
	 label dimList = 1;
     min_.setSize(entries_.size());
     dim_.setSize(entries_.size());
	 delta_.setSize(entries_.size());
     max_.setSize(entries_.size());
		
	 label index = 0;
     forAll(entries_,i)
     {
        dim_[i] = readLabel(entries_[i].lookup("N"));
	    max_[i] = readScalar(entries_[i].lookup("max"));
    	min_[i] = readScalar(entries_[i].lookup("min"));
		delta_[i] = (max_[i] - min_[i]) / dim_[i];
        dimList *= dim_[i];
		fieldIndices_.insert(entries_[i].lookup("name"),index);
		index ++;
    }
	
	forAll(output_,i)
    {
		fieldIndices_.insert(output_[i].lookup("name"),index);
		index++;
	}	
		
    List<List<Type> >& internal = *this;
				
	internal.setSize(dimList);
		
	forAll(internal, i)
	{
		internal[i].setSize(entries_.size()+output_.size());
	}	
}

template<class Type>
void Foam::interpolationLookUpTable<Type>::readTable()
{
    // preserve the original (unexpanded) fileName to avoid absolute paths
    // appearing subsequently in the write() method
    
	fileName fName("constant"/fileName_);
	
    fName.expand();
	
	IFstream  controlFile(fName);
	
	const dictionary control(controlFile);
	
    control.lookup("fields") >> entries_;
	
	control.lookup("output") >> output_;

	control.lookup("values") >> *this;
	
	dimensionTable();
	
    check();

    if (this->size() == 0)
    {
        FatalErrorIn
        (
                "Foam::interpolationLookUpTable<Type>::readTable()"
        )   << "table is empty" << nl
            << exit(FatalError);
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable()
:
	List<List<Type> >(),
    boundsHandling_(interpolationLookUpTable::CLAMP),
    fileName_("fileNameIsUndefined")
{}


template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable
(
	const fileName& fn
)
:
	List<List<Type> >(),
    boundsHandling_(interpolationLookUpTable::CLAMP),
    fileName_(fn),
	dim_(0),
	min_(0),
	delta_(0.),
	max_(0.),
	entries_(0),
	output_(0)
{
    readTable();
}

template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable
(
     const interpolationLookUpTable& interpTable
)
:
	List<List<Type> >(interpTable),
    boundsHandling_(interpTable.boundsHandling_),
    fileName_(interpTable.fileName_)
{}

template<class Type>
Foam::interpolationLookUpTable<Type>::interpolationLookUpTable
(
	const dictionary& dict
)
:
	List<List<Type> >(),
    boundsHandling_(interpolationLookUpTable::CLAMP),
    fileName_(fileName(dict.lookup("fileName")).expand()),
    dim_(0),
    min_(.0),
    delta_(.0),
    max_(.0),
    entries_(dict.lookup("fields")),
	output_(dict.lookup("output")),
	fieldIndices_(0)
    {
		dimensionTable();
    }

		
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::interpolationLookUpTable<Type>::insertFields
(
	List<dictionary>& fields,
 	List<dictionary>& output
) 
{
	entries_ = fields;
	output_  = output;
	dimensionTable();
}
template<class Type>
Foam::word Foam::interpolationLookUpTable<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case interpolationLookUpTable::ERROR:
        {
            enumName = "error";
            break;
        }
        case interpolationLookUpTable::WARN:
        {
            enumName = "warn";
            break;
        }
        case interpolationLookUpTable::CLAMP:
        {
            enumName = "clamp";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::interpolationLookUpTable<Type>::boundsHandling
Foam::interpolationLookUpTable<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return interpolationLookUpTable::ERROR;
    }
    else if (bound == "warn")
    {
        return interpolationLookUpTable::WARN;
    }
    else if (bound == "clamp")
    {
        return interpolationLookUpTable::CLAMP;
    }
    else
    {
        WarningIn
        (
                "Foam::interpolationLookUpTable<Type>::wordToBoundsHandling(const word&)"
        )   << "bad outOfBounds specifier " << bound << " using 'warn'" << endl;

        return interpolationLookUpTable::WARN;
    }
}


template<class Type>
typename Foam::interpolationLookUpTable<Type>::boundsHandling
Foam::interpolationLookUpTable<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::interpolationLookUpTable<Type>::check() const
{
//  check order in the first dimension.
	
    scalar prevValue = List<List<Type> >::operator[](0).operator[](0);
    label dim = 1 ;
    for(int j = 1 ;j <=  dim_.size()-1 ; j++) dim *=dim_[j];

    for (label i = 1; i< dim_[0]; ++i)
    {
        label index = i*dim + 1;
        const scalar currValue =
        List<List<Type> >::operator[](index).operator[](0);
        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorIn
            (
                    "Foam::interpolationLookUpTable<Type>::checkOrder() const"
            )   << "out-of-order value: "
                << currValue << " at index " << index << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}

template<class Type>
void Foam::interpolationLookUpTable<Type>::write(Ostream& os) const
{
   os.writeKeyword("fields");
   os << entries_ << token::END_STATEMENT << nl;

   os.writeKeyword("output");
   os << output_ << token::END_STATEMENT << nl;
   
    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Foam::interpolationTable<Type>::write()"
        )   << "table is empty" << nl
            << exit(FatalError);
    }
	os.writeKeyword("values");
    os << *this << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type> 
Foam::List<Type>&  
Foam::interpolationLookUpTable<Type>::operator[]
(
	const label i
)
{
    label ii = i;
    label n  = this->size();

    if (n <= 1)
    {
        ii = 0;
    }
    else if (ii < 0)
    {
        switch (boundsHandling_)
        {
            case interpolationLookUpTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolationLookUpTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationLookUpTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolationLookUpTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolationLookUpTable::CLAMP:
            {
                ii = 0;
                break;
            }
        }
    }
    else if (ii >= n)
    {
        switch (boundsHandling_)
        {
            case interpolationLookUpTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolationLookUpTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationLookUpTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolationLookUpTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolationLookUpTable::CLAMP:
            {
                ii = n - 1;
                break;
            }
        }
    }

	return List<List<scalar> >::operator[](ii);
	       
}

template<class Type>
const Foam::List<Type>& 
Foam::interpolationLookUpTable<Type>::operator[]
(
	const label i
) const
{
//    label index = findField(field);
	label ii = i;
    label n  = this->size();

    if (n <= 1)
    {
        ii = 0;
    }
    else if (ii < 0)
    {
        switch (boundsHandling_)
        {
            case interpolationLookUpTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolationLookUpTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationLookUpTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolationLookUpTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolationLookUpTable::CLAMP:
            {
                ii = 0;
                break;
            }
        }
    }
    else if (ii >= n)
    {
        switch (boundsHandling_)
        {
            case interpolationLookUpTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolationLookUpTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationLookUpTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolationLookUpTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolationLookUpTable::CLAMP:
            {
                ii = n - 1;
                break;
            }
        }
    }

	return List<List<scalar> >::operator[](ii);
}

template<class Type> 
Foam::scalar
Foam::interpolationLookUpTable<Type>::operator()
(
	scalarList& ivector, const word& ifield
)const
{
	label n = this->size();
	
    if (n <= 1)
    {
		FatalErrorIn
        (
         	"Foam::interpolationLookUpTable<Type>::operator[]"
             "(scalarList& vector, const word& ifield) const"
        )   << "empty table" << nl
        	<< exit(FatalError);
        
    }

	forAll(ivector,i)
    {
        scalar lookupValue = ivector[i]; 
        scalar minLimit = min_[i];
        scalar maxLimit = max_[i];
		
        if (lookupValue < minLimit)
        {
            switch (boundsHandling_)
            {
                case interpolationLookUpTable::ERROR:
                {
                    FatalErrorIn
                    (
                        "Foam::interpolationLookUpTable<Type>::operator[]"
                        "(const scalar) const"
                    )   << "value (" << lookupValue << ") underflow" << nl
                        << exit(FatalError);
                    break;
                }
                case interpolationLookUpTable::WARN:
                {
                    WarningIn
                    (
                        "Foam::interpolationLookUpTable<Type>::operator[]"
                        "(const scalar) const"
                    )   << "value (" << lookupValue << ") underflow" << nl
                        << "    Continuing with the first entry"
                        << endl;
                    // fall-through to 'CLAMP'
                }
                case interpolationLookUpTable::CLAMP:
                {
				 
					ivector[i] = min_[i];
//               return List<FixedList<Type,3> >::operator[](0).operator[](i);
                    break;
                }
            }
        }
        else if (lookupValue >= maxLimit)
        {
            switch (boundsHandling_)
            {
                case interpolationLookUpTable::ERROR:
                {
                    FatalErrorIn
                    (
                    "Foam::interpolationLookUpTable<Type>::operator[]"
                    "(const label) const"
                    )   << "value (" << lookupValue << ") overflow" << nl
                        << exit(FatalError);
                    break;
                }
                case interpolationLookUpTable::WARN:
                {
                    WarningIn
                    (
                        "Foam::interpolationLookUpTable<Type>::operator[]"
                        "(const label) const"
                    )   << "value (" << lookupValue << ") overflow" << nl
                        << "    Continuing with the last entry"
                        << endl;
                // fall-through to 'CLAMP'
                }
                case interpolationLookUpTable::CLAMP:
                {	
				
					ivector[i] = max_[i];
//                    return List<FixedList<Type ,3> >::operator[](n-1).operator[](i);
                    break;
                }
            }
        }				
	}

	label fieldIndex = findField(ifield);	
	
	label lo = index(ivector);
	
    label dim = 1;
    for(int j = 1 ;j <  dim_.size() ; j++)
	{
		dim *=dim_[j];
	}
    
	label hi = lo + dim;
	
//	Info << "lo:" << lo << endl;
//	Info << "hi:" << hi << endl;
	
//	dimension to interpolate.
	
	scalar lookupValue = ivector[0];
	
    List<scalar> outputList
    (
		List<List<Type> >::operator[](lo)
     +  (
        	List<List<Type> >::operator[](hi)
       	  - List<List<Type> >::operator[](lo)
        )
      * (
        	lookupValue
          - List<List<Type> >::operator[](lo).operator[](0)
        )
       /(
        	List<List<Type> >::operator[](hi).operator[](0)
          - List<List<Type> >::operator[](lo).operator[](0)
        )
     );
	return(outputList[fieldIndex]);
}

template<class Type>
inline Foam::scalar
Foam::interpolationLookUpTable<Type>::operator()
(
	const scalar T,
	const label ifield
)const
{
	scalarList retvals(1,T);
	
	label lo = index(retvals);
	
	label dim = 1;
	for(int j = 1 ;j <  dim_.size() ; j++)
	{
		dim *=dim_[j];
	}
    
	label hi = lo + dim;
	
	scalar lookupValue = T;
	
	scalar outputScalar
			(
			List<List<Type> >::operator[](lo).operator[](ifield)
			+  (
			List<List<Type> >::operator[](hi).operator[](ifield)
			- List<List<Type> >::operator[](lo).operator[](ifield)
			   )
			* (
			lookupValue
			- List<List<Type> >::operator[](lo).operator[](0)
			  )
       /(
			List<List<Type> >::operator[](hi).operator[](0)
			- List<List<Type> >::operator[](lo).operator[](0)
		)
			);
//	return(List<List<Type> >::operator[](lo).operator[](ifield));
	return(outputScalar);
}


template<class Type>
inline Foam::scalar
Foam::interpolationLookUpTable<Type>::operator()
(
const scalar ft,
const scalar ftVar,
const label ifield
)const
{
	scalarList retvals(2,0.);
	retvals[0] = ft;
	retvals[1] = ftVar;
	
	label lo = index(retvals);
	
	label dim = 1;
	for(int j = 1 ;j <  dim_.size() ; j++)
	{
		dim *=dim_[j];
	}
    
	label hi = lo + dim;
	
	scalar lookupValue = ft;
	
	scalar outputScalar
			(
			List<List<Type> >::operator[](lo).operator[](ifield)
			+  (
			List<List<Type> >::operator[](hi).operator[](ifield)
			- List<List<Type> >::operator[](lo).operator[](ifield)
			   )
			* (
			lookupValue
			- List<List<Type> >::operator[](lo).operator[](0)
			  )
       /(
			List<List<Type> >::operator[](hi).operator[](0)
			- List<List<Type> >::operator[](lo).operator[](0)
		)
			);
	return(outputScalar);
}


template<class Type>
Foam::interpolationLookUpTable<Type>&
Foam::interpolationLookUpTable<Type>::operator*=
(
 	const scalar s
)
{
	List<List<Type> >& internal = *this;
	forAll(output_,i)
    {
		label index = findField(output_[i].lookup("name"));
		forAll(internal, i)
		{
			internal[i][index] *= s ;
		}
	}
	
	return *this;
}


template<class Type>
Foam::interpolationLookUpTable<Type>&
Foam::interpolationLookUpTable<Type>::operator=
(
	const interpolationLookUpTable<Type>& itp
)
{
	List<List<Type> >& internal = *this;
	
	forAll(output_,i)
    {
		label index = findField(output_[i].lookup("name"));
		
		forAll(internal, i)
		{
			internal[i][index] = itp[i][index];
		}
	}
	return *this;
}

template<class Type>
Foam::interpolationLookUpTable<Type>&
Foam::interpolationLookUpTable<Type>::operator+= //
(
	const interpolationLookUpTable<Type>& itp1
)
{
	List<List<Type> >& internal = *this;
	forAll(output_,i)
    {
		label index = findField(output_[i].lookup("name"));

		forAll(itp1, i)
		{
			internal[i][index] += itp1[i][index];
		}
	}
	return *this;
}

template<class Type>
Foam::interpolationLookUpTable<Type>&
Foam::interpolationLookUpTable<Type>::operator-= //
(
	const interpolationLookUpTable<Type>& itp1
)
{
	List<List<Type> >& internal = *this;
	forAll(output_,i)
    {
		label index = findField(output_[i].lookup("name"));

		forAll(itp1, i)
		{
			internal[i][index] -= itp1[i][index];
		}
	}
	return *this;
}




