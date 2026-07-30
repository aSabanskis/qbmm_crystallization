/*---------------------------------------------------------------------------*\
  Solution crystallization nucleation model for OpenQBMM.
\*---------------------------------------------------------------------------*/

#ifndef solutionNucleation_H
#define solutionNucleation_H

#include "nucleationModel.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

class solutionNucleation
:
    public nucleationModel
{
private:

    // Empirical primary-nucleation coefficient in the Liu model's unit system
    scalar kS_;

    // Activation energy and gas constant
    dimensionedScalar Es_;
    dimensionedScalar R_;

    // Size at which new crystals enter the population balance
    scalar nucleationSize_;

    // Concentration and temperature fields
    const volScalarField& C_;
    const volScalarField& T_;

    scalar saturationConcentration(const scalar T) const;

public:

    TypeName("solution");

    solutionNucleation
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~solutionNucleation();

    virtual scalar nucleationSource
    (
        const label& momentOrder,
        const label celli,
        const label environment = 0
    ) const;
};

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

#endif

// ************************************************************************* //
