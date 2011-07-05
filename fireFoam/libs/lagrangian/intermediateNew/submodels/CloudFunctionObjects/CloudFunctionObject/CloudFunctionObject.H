/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::CloudFunctionObject

Description
    Templated cloud functions

SourceFiles
    CloudFunctionObject.C
    NewCloudFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef CloudFunctionObject_H
#define CloudFunctionObject_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class CloudFunctionObject Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudFunctionObject
{
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


protected:

    // Protected member functions

        //- Write post-processing info
        virtual void write() = 0;


public:

    //- Runtime type information
    TypeName("CloudFunctionObject");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CloudFunctionObject,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        CloudFunctionObject(CloudType& owner);

        //- Construct from dictionary
        CloudFunctionObject
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~CloudFunctionObject();


    //- Selectors
    static autoPtr<CloudFunctionObject<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );
    static autoPtr<CloudFunctionObject<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner,
        const word& modelType
    );


    // Member Functions

        // Access

            //- Return the owner cloud dictionary
            inline const dictionary& dict() const;

            //- Return const access the owner cloud object
            inline const CloudType& owner() const;

            //- Return non-const access the owner cloud object for manipulation
            inline CloudType& owner();

            //- Return the coefficients dictionary
            inline const dictionary& coeffDict() const;


        // Evaluation

            //- Post-evolve hook
            virtual void postEvolve();

            //- Post-patch hook
            virtual void postPatch
            (
                const typename CloudType::parcelType& p,
                const label patchI
            ) = 0;

            //- Post-face hook
            virtual void postFace(const typename CloudType::parcelType& p) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCloudFunctionObject(CloudType)                                    \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(CloudFunctionObject<CloudType>, 0);   \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        CloudFunctionObject<CloudType>,                                       \
        dictionary                                                            \
    );


#define makeCloudFunctionObjectType(SS, CloudType, ParcelType)                \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    CloudFunctionObject<CloudType<ParcelType> >::                             \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


#define makeCloudFunctionObjectThermoType(SS, CloudType, ParcelType, ThermoType)\
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    CloudFunctionObject<CloudType<ParcelType<ThermoType> > >::                \
        adddictionaryConstructorToTable                                       \
            <SS<CloudType<ParcelType<ThermoType> > > >                        \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CloudFunctionObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CloudFunctionObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
