IOdictionary crystallizationProperties
(
    IOobject
    (
        "crystallizationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar DC("DC", dimArea/dimTime, crystallizationProperties);
const scalar ScT = crystallizationProperties.lookupOrDefault<scalar>("ScT", 0.7);

Info<< "\nRight solver!\n" << endl;
dimensionedScalar rhoCrystal
(
    "rhoCrystal",
    dimMass/dimVolume,
    crystallizationProperties
);

const scalar kv =
    crystallizationProperties.lookupOrDefault<scalar>("kv", 1.0);

const bool limitC =
    crystallizationProperties.lookupOrDefault<bool>("limitC", true);


volScalarField DCEff
(
    IOobject
    (
        "DCEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DC + turbulence->nut()/ScT
);

fvScalarMatrix CEqn
(
    fvm::ddt(rho, C)
  + fvm::div(phi, C)
  - fvm::laplacian(rho*DCEff, C)
);

CEqn.relax();
CEqn.solve();

C.correctBoundaryConditions();


const volScalarField& m3 =
    mesh.lookupObject<volScalarField>
    (
        "moment.3.populationBalance"
    );

volScalarField m3Old("m3Old", m3);

populationBalance->solve();

volScalarField dm3("dm3", m3 - m3Old);

C.primitiveFieldRef() -= rhoCrystal.value()*kv*dm3.primitiveField();

C.primitiveFieldRef() -= rhoCrystal.value()*kv*dm3.primitiveField();

if (limitC)
{
    C.max(0.0);
}

C.correctBoundaryConditions();
