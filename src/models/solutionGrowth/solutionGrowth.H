/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015 by Matteo Icardi and 2017 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
-------------------------------------------------------------------------------
2017-03-28 Alberto Passalacqua: Adapted to single scalar calculation.
2020-03-12 Alberto Passalacqua: Added access to fvMesh.
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::populationBalanceSubModels::growthModels::solutionGrowth

Description
    Growth from solution

SourceFiles
    solutionGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef solutionGrowth_H
#define solutionGrowth_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                    Class solutionGrowth Declaration
\*---------------------------------------------------------------------------*/

class solutionGrowth
:
    public growthModel
{
private:

    // Private data

        //- Minimum value of the abscissa to apply the model
        scalar minAbscissa_;

        //- Maximum value of the abscissa to apply the model
        scalar maxAbscissa_;

        dimensionedScalar kG_;
        dimensionedScalar Csat0_;
        dimensionedScalar dCsatdT_;
        dimensionedScalar Tref_;
        scalar g_;
        dimensionedScalar Eg_;
        dimensionedScalar R_;
        dimensionedScalar alphaG_;
        scalar Lg_;

public:

        //- Runtime type information
        TypeName("solution");


    // Constructors

        //- Construct from components
        solutionGrowth
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~solutionGrowth();


    scalar saturationConcentration
    (
        const scalar T
    ) const;

    scalar growthRate
    (
        const scalar L,
        const scalar C,
        const scalar T
    ) const;
    
    // Member Functions

        //- Growth model
        virtual scalar Kg
        (
            const scalar& abscissa,
            const bool lengthBased,
            const label environment = 0
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
