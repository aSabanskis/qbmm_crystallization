Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

Info<< "Reading concentration C\n" << endl;

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

pressureControl pressureControl(p, rho, pimple.dict(), false);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

mesh.setFluxRequired(p_rgh.name());

label pRefCell = 0;
scalar pRefValue = 0.0;

if (p_rgh.needReference())
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

dimensionedScalar initialMass("initialMass", fvc::domainIntegrate(rho));

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"
#include "createRadiationModel.H"
#include "createFvOptions.H"


const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

Info<< "Reading crystallizationProperties\n" << endl;

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

dimensionedScalar rhoCrystal
(
    "rhoCrystal",
    dimMass/dimVolume,
    crystallizationProperties
);

// Solubility polynomial, Csat(Tc) with Tc in degC, highest order first.  It has
// to match the one the growth and nucleation models use, so that the solute
// balance and the kinetics agree on where equilibrium is.
const scalarList CsatCoeffs
(
    crystallizationProperties.get<scalarList>("CsatCoeffs")
);

if (CsatCoeffs.size() != 5)
{
    FatalIOErrorInFunction(crystallizationProperties)
        << "CsatCoeffs must be a quartic in T [degC], i.e. 5 coefficients, "
        << "but " << CsatCoeffs.size() << " were given"
        << exit(FatalIOError);
}

auto saturation = [](const scalar T, const scalarList& a) -> scalar
{
    const scalar Tc = T - 273.15;
    return ((((a[0]*Tc + a[1])*Tc + a[2])*Tc + a[3])*Tc + a[4]);
};

Info<< "Reading populationBalanceProperties\n" << endl;

IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

surfaceScalarField phiByRho("phiByRho", phi/fvc::interpolate(rho));

autoPtr<populationBalanceModel> populationBalance
(
    populationBalanceModel::New
    (
        "populationBalance", populationBalanceProperties, phiByRho
    )
);
