/*---------------------------------------------------------------------------*\
  Solution crystallization nucleation model for OpenQBMM.
\*---------------------------------------------------------------------------*/

#include "solutionNucleation.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(solutionNucleation, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        solutionNucleation,
        dictionary
    );
}
}
}


Foam::populationBalanceSubModels::nucleationModels::solutionNucleation
::solutionNucleation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    kS_(dict.get<scalar>("kS")),
    Es_("Es", dimEnergy/dimMoles, dict),
    R_("R", dimEnergy/dimMoles/dimTemperature, dict),
    nucleationSize_(dict.lookupOrDefault<scalar>("nucleationSize", 0)),
    C_(mesh.lookupObject<volScalarField>(dict.lookupOrDefault<word>("C", "C"))),
    T_(mesh.lookupObject<volScalarField>(dict.lookupOrDefault<word>("T", "T")))
{}


Foam::populationBalanceSubModels::nucleationModels::solutionNucleation
::~solutionNucleation()
{}


Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::solutionNucleation
::saturationConcentration(const scalar T) const
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
Foam::populationBalanceSubModels::nucleationModels::solutionNucleation
::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const
{
    const scalar supersaturation =
        max(C_[celli] - saturationConcentration(T_[celli]), scalar(0));

    const scalar B =
        kS_
       *supersaturation
       *Foam::exp(-Es_.value()/(R_.value()*T_[celli]));

    return B*Foam::pow(nucleationSize_, momentOrder);
}

// ************************************************************************* //
