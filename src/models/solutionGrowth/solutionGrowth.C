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

\*---------------------------------------------------------------------------*/

#include "solutionGrowth.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(solutionGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        solutionGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::solutionGrowth
::solutionGrowth
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    growthModel(dict, mesh),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", scalar(0))),
    maxAbscissa_(dict.lookupOrDefault("maxAbscissa", GREAT)),
    kG_("kG", dimLength/dimTime, dict),
    Csat0_("Csat0", dimMass/dimVolume, dict),
    dCsatdT_("dCsatdT", dimMass/dimVolume/dimTemperature, dict),
    Tref_("Tref", dimTemperature, dict),
    g_(dict.lookupOrDefault("g", scalar(1))),
    Eg_("Eg", dimEnergy/dimMoles, dict),
    R_("R", dimEnergy/dimMoles/dimTemperature, dict),
    alphaG_("alphaG", inv(dimLength), dict),
    Lg_(dict.lookupOrDefault("Lg", scalar(1)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::solutionGrowth
::~solutionGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::growthModels::solutionGrowth
::saturationConcentration
(
    const scalar T
) const
{
    const scalar TC = T - scalar(273.15);

    return
          4e-5*Foam::pow(TC, 4)
        - 0.0034*Foam::pow(TC, 3)
        + 0.1024*Foam::pow(TC, 2)
        - 0.6255*TC
        + 13.237;
}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::solutionGrowth
::growthRate
(
    const scalar L,
    const scalar C,
    const scalar T
) const
{
    const scalar Csat = saturationConcentration(T);

    const scalar deltaC =
        max(C - Csat, scalar(0));
    return
        kG_.value()
       *Foam::pow(deltaC, g_)
       *Foam::exp(-Eg_.value()/(R_.value()*T))
       *Foam::pow(1.0 + alphaG_.value()*L, Lg_);
}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::solutionGrowth::Kg
(
    const scalar& abscissa,
    const bool lengthBased,
    const label environment
) const
{
    const volScalarField& C =
        mesh_.lookupObject<volScalarField>("C");

    const volScalarField& T =
        mesh_.lookupObject<volScalarField>("T");

    return
        growthRate
        (
            abscissa,
            C[environment],
            T[environment]
        )
       *pos0(abscissa - minAbscissa_)
       *neg0(abscissa - maxAbscissa_);
}

// ************************************************************************* //
