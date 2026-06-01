/*---------------------------------------------------------------------------*\
License
    This file is part of interIsoThermalFoam.

    interIsoThermalFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    interIsoThermalFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with interIsoThermalFoam. 


Authors
    Zhao Wang
\*---------------------------------------------------------------------------*/
#include "interIsoThermalEnergy.H"

namespace Foam
{

interIsoThermalEnergy::interIsoThermalEnergy
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const volScalarField& rho,
    const surfaceScalarField& rhoPhi,
    const dimensionedScalar& rho1,
    const dimensionedScalar& rho2,
    fv::options& fvOptions
)
:
    // Link external objects to members in this class
    mesh_(mesh),
    alpha1_(alpha1),
    alpha2_(alpha2),
    rho_(rho),
    rhoPhi_(rhoPhi),
    rho1_(rho1),
    rho2_(rho2),
    fvOptions_(fvOptions),

    thermalProperties_
    (
        IOobject
        (
            "thermalProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    cp1_
    (
        "cp",
        dimensionSet(0, 2, -2, -1, 0, 0, 0),
        thermalProperties_.subDict("phase1")
    ),

    cp2_
    (
        "cp",
        dimensionSet(0, 2, -2, -1, 0, 0, 0),
        thermalProperties_.subDict("phase2")
    ),

    kappa1_
    (
        "kappa",
        dimensionSet(1, 1, -3, -1, 0, 0, 0),
        thermalProperties_.subDict("phase1")
    ),

    kappa2_
    (
        "kappa",
        dimensionSet(1, 1, -3, -1, 0, 0, 0),
        thermalProperties_.subDict("phase2")
    ),

    Prt_
    (
        "Prt",
        dimless,
        thermalProperties_.lookupOrDefault<scalar>("Prt", 0.85)
    ),

    Tsat_
    (
        "Tsat",
        dimensionSet(0, 0, 0, 1, 0, 0, 0),
        thermalProperties_.subDict("phaseChange")
    ),

    Lv_
    (
        "Lv",
        dimensionSet(0, 2, -2, 0, 0, 0, 0),
        thermalProperties_.subDict("phaseChange")
    ),

    Cevap_
    (
        "Cevap",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        thermalProperties_.subDict("phaseChange")
    ),

    Ccond_
    (
        "Ccond",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        thermalProperties_.subDict("phaseChange")
    ),

    T_
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    cp_ // Specific heat of mixture
    (
        IOobject
        (
            "cp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "cp",
            dimensionSet(0, 2, -2, -1, 0, 0, 0),
            0.0
        )
    ),

    kappa_
    (
        IOobject
        (
            "kappa",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "kappa",
            dimensionSet(1, 1, -3, -1, 0, 0, 0),
            0.0
        )
    ),

    rhoCp_
    (
        IOobject
        (
            "rhoCp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "rhoCp",
            dimensionSet(1, -1, -2, -1, 0, 0, 0),
            0.0
        )
    ),

    rhophicp_
    (
        IOobject
        (
            "rhophicp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "rhophicp",
            dimensionSet(1, 2, -3, -1, 0, 0, 0),
            0.0
        )
    ),

    mDotEvap_
    (
        IOobject
        (
            "mDotEvap", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("mDotEvap", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),

    mDotCond_
    (
        IOobject
        (
            "mDotCond", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("mDotCond", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),

    mDot_
    (
        IOobject
        (
            "mDot", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("mDot", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),

    latentSource_
    (
        IOobject
        (
            "latentSource", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("latentSource", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)
    ),

    vDot_
    (
        IOobject
        (
            "vDot", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("vDot", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),

    SuAlpha_
    (
        IOobject
        (
            "SuAlpha", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("SuAlpha", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),

    SpAlpha_
    (
        IOobject
        (
            "SpAlpha", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("SpAlpha", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    )

{
    correct();
    rhoCp_.oldTime();
}

void interIsoThermalEnergy::correct()
{
    rhoCp_ = alpha1_*rho1_*cp1_ + alpha2_*rho2_*cp2_;
    rhoCp_.correctBoundaryConditions();

    cp_ = rhoCp_/rho_;
    cp_.correctBoundaryConditions();

    kappa_ = alpha1_*kappa1_ + alpha2_*kappa2_;
    kappa_.correctBoundaryConditions();

    
    rhophicp_ = fvc::interpolate(cp_)*rhoPhi_;
    rhophicp_.correctBoundaryConditions();

    correctPhaseChange();
}

void interIsoThermalEnergy::correctPhaseChange() // Lee's phase-change model
{
    const dimensionedScalar zero("zero", dimless, 0.0);

    tmp<volScalarField> tSuperHeat = max((T_ - Tsat_)/Tsat_, zero);
    tmp<volScalarField> tSubCool = max((Tsat_ - T_)/Tsat_, zero);

    const volScalarField& superHeat = tSuperHeat();
    const volScalarField& subCool = tSubCool();

    mDotEvap_ = Cevap_*rho1_*alpha1_*superHeat;
    mDotCond_ = Ccond_*rho2_*alpha2_*subCool;
    mDot_ = mDotEvap_ - mDotCond_;
    latentSource_ = -mDot_*Lv_; // Temperature eqn source
    vDot_ = mDot_*(1.0/rho2_ - 1.0/rho1_); //Pressure eqn source

    SuAlpha_ = Ccond_*(rho2_/rho1_)*subCool;
    SpAlpha_ = - ( Cevap_*superHeat + Ccond_*(rho2_/rho1_)*subCool);

    mDotEvap_.correctBoundaryConditions();
    mDotCond_.correctBoundaryConditions();
    mDot_.correctBoundaryConditions();
    latentSource_.correctBoundaryConditions();
    vDot_.correctBoundaryConditions();
    SuAlpha_.correctBoundaryConditions();
    SpAlpha_.correctBoundaryConditions();
}

void interIsoThermalEnergy::solve()
{
    // Update properties and source terms
    correct();

    volScalarField kappaEff
    (
        IOobject
        (
            "kappaEff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        kappa_
    );

    if (mesh_.foundObject<volScalarField>("nut"))
    {
        const volScalarField& nut = mesh_.lookupObject<volScalarField>("nut");

        kappaEff += rhoCp_*nut/Prt_;
    }

    volScalarField latentCoeff
    (
        IOobject
        (
            "latentCoeff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (
            Ccond_*rho2_*alpha2_*pos(Tsat_ - T_)
          + Cevap_*rho1_*alpha1_*pos(T_ - Tsat_)
        )*Lv_/Tsat_
    );

    fvScalarMatrix TEqn
    (
  
        fvm::ddt(rhoCp_, T_)

      + fvm::div(rhophicp_, T_)

      - fvm::Sp(fvc::ddt(rhoCp_) + fvc::div(rhophicp_), T_)

      - fvm::laplacian(kappaEff, T_)

     ==
      - fvm::Sp(latentCoeff, T_)

      + latentCoeff*Tsat_

      + fvOptions_(rhoCp_, T_)
    );

    TEqn.relax();

    fvOptions_.constrain(TEqn);

    TEqn.solve();

    fvOptions_.correct(T_);

    T_.correctBoundaryConditions();

    correctPhaseChange();
}

const volScalarField& interIsoThermalEnergy::T() const
{
    return T_;
}


const volScalarField& interIsoThermalEnergy::mDot() const
{
    return mDot_;
}


const volScalarField& interIsoThermalEnergy::latentSource() const
{
    return latentSource_;
}


const volScalarField& interIsoThermalEnergy::vDot() const
{
    return vDot_;
}


const volScalarField& interIsoThermalEnergy::SuAlpha() const
{
    return SuAlpha_;
}


const volScalarField& interIsoThermalEnergy::SpAlpha() const
{
    return SpAlpha_;
}

}



