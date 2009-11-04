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

#include "fireMixture.H"
#include "fvMesh.H"

#include "IFstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ThermoType>
const char* fireMixture<ThermoType>::specieNames_[6] = 
{
    "ft", 
    "b",
    "ftVar",
    "fu",
    "ox",
    "pr"
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::fireMixture<ThermoType>::fireMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
:
    basicMultiComponentMixture
    (
        thermoDict, 
        speciesTable(nSpecies_, specieNames_), 
        mesh
    ),

    stoicRatio_(thermoDict.lookup("stoichiometricAirFuelMassRatio")),
    fuel_(thermoDict.lookup("fuel")),
    oxidant_(thermoDict.lookup("oxidant")),
    products_(thermoDict.lookup("burntProducts")),
    mixture_("mixture", fuel_),
    ft_(Y("ft")),
    b_(Y("b")),
    ftVar_(Y("ftVar")),
    fu_(Y("fu")),
    ox_(Y("ox")),
    pr_(Y("pr")),
    //LookUpTable_(),
    pdfMethod_()
{
    const dictionary dictBetaIntegration =
    ft_.mesh().solutionDict().subDict("betaIntegration");

    pdfMethod_ = readLabel(dictBetaIntegration.lookup("pdfMethod"));

    Info << "pdfMethod = " << pdfMethod_ << endl;

    if (pdfMethod_ == fireMixture::lookUpTable)
    {
        //LookUpTable_ = fileName(dictBetaIntegration.lookup("LookUpTable"));
    }

    if(pdfMethod_ == fireMixture::recursive)
    {
        //Read coeff_[]
        fileName fName(ft_.time().constant()/"betaPdfCoeff.txt"); 

        fName.expand();

        IFstream betaPdfCoeffFile(fName); 

        if (betaPdfCoeffFile.good())
        {
            for (int i=0; i<nPoly_; i++)
            {
                betaPdfCoeffFile >> coeff_[i];
                Info << "betaPdfCoeff[" << i << "] = " << coeff_[i] << endl;
            }
        }
        else
        {
            FatalErrorIn("fireMxiture")
                << "Cannot read file" << fName
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& fireMixture<ThermoType>::mixture
(
    const scalar ft,
    const scalar b,
    const scalar ftVar,
    scalar& fu,
    scalar& ox,
    scalar& pr
) const
{

    if (ft < 0.0001)     //for nonreacting zone, save the time of doing beta pdf
    {
        fu = 0.0;
        ox = 1.0;
        pr = 0.0;
        return oxidant_; //if ft is unbounded, do clip
    }
    else if(ft > 0.9999)   //for nonreacting zone, save the time of doing beta pdf
    {
        fu = 1.0;
        ox = 0.0;
        pr = 0.0;
        return fuel_;    //if ft is unbounded, do clip
    }
    else
    {
        if (pdfMethod_ == fireMixture::lookUpTable)
        {
            FatalErrorIn("fireMixture::mixture() const")
            << " lookUpTable option was turned off to be compatible with FVM radiation"
            << exit(FatalError);

            //fu = LookUpTable_(ft, ftVar, 2);
            //fu = b*ft + (1.0 - b)*fu;
        }

        if (pdfMethod_ == fireMixture::recursive)
        {
            fu = getMassFraction(ft, ftVar, 1.0 - b);
        }

        if (pdfMethod_ == fireMixture::delta)
        {
            fu = b*ft + (1.0 - b)*fres(ft, stoicRatio().value());
        }

        ox = 1 - ft - (ft - fu)*stoicRatio().value();
        pr = 1 - fu - ox;

        mixture_ = fu/fuel_.W()*fuel_;
        mixture_ += ox/oxidant_.W()*oxidant_;
        mixture_ += pr/products_.W()*products_;

        return mixture_;
    }
}

template<class ThermoType>
void fireMixture<ThermoType>::read(const dictionary& thermoDict)
{
    stoicRatio_ = thermoDict.lookup("stoichiometricAirFuelMassRatio");
    fuel_ = ThermoType(thermoDict.lookup("fuel"));
    oxidant_ = ThermoType(thermoDict.lookup("oxidant"));
    products_ = ThermoType(thermoDict.lookup("burntProducts"));
}


template<class ThermoType>
scalar fireMixture<ThermoType>::getMassFraction
(
    const scalar mean, //filtered mixture fraction
    const scalar var,  //SGS variance of mixture fraction
    const scalar c     //progress variable (c=1-b)
) const
{
    scalar fu;         //fuel mass fraction
    scalar ft;         //mixture fraction
    scalar fu_m;

//    const scalar epsilon = 1e-10;
    const scalar epsilon = 1e-5;

    //  first check if var is within the theoretical limit [0,mean*(1-mean)] 
    if (c < 1e-10)            // no reaction, mixing solution 
    {
        fu_m = mean;
    }

    else if (var <= epsilon)  // use Dirac Delta
    {
        ft = mean;
        fu = (1.0-c)*ft + c*fres(ft, stoicRatio().value());
        fu_m = fu;
    }

    else if ( var >= mean*(1.0-mean) - epsilon)  // use Bernoulli distribution
    {
        fu_m = mean;
    }

    else                //use beta distribution
    {
    	scalar a = mean*(mean*(1.0-mean)/var-1.0);   
    	scalar b = (1.0-mean)/mean*a;               

        if ( a < 0 or b < 0 ) 
        {
          Info << "mean = " << mean << endl;
          Info << "var = " << var << endl;
          FatalErrorIn("fireMxiture::getMassFraction_recu")
              << "Error: in Beta PDF, a = " << a << "b = " << b 
              << abort(FatalError);
	}

        fu_m = coeff_[0];
	scalar factor = 1.0;
        for (int i=1; i<nPoly_; i++)
        {
          factor *= (a+scalar(i-1))/(a+b+scalar(i-1));
          fu_m += coeff_[i]*factor;
	}

        //include c effect
        fu_m = (1.0 - c)*mean + c*fu_m;

        //check
        if (fu_m < 0.0)
        {
          fu_m = 0.0;
        }
        else if (fu_m > 1.0) 
        {
          fu_m = 1.0;
        }

        //check if fu_m >= fu
        ft = mean;
        fu = (1.0 - c)*ft + c*fres(ft, stoicRatio().value());
        if (fu_m < fu)
        {
          fu_m = fu;
	}
        else if (fu_m > mean)
        {
          fu_m = mean;
        } 
    }
    return fu_m;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
