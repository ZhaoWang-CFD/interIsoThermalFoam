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

Class
    interIsoThermalEnergy


Authors
    Zhao Wang
\*---------------------------------------------------------------------------*/

#ifndef interIsoThermalEnergy_H
#define interIsoThermalEnergy_H

#include "fvCFD.H"
#include "fvOptions.H"

namespace Foam
{

class interIsoThermalEnergy
{
    const fvMesh& mesh_;

    const volScalarField& alpha1_;
    const volScalarField& alpha2_;
    const volScalarField& rho_;
    const surfaceScalarField& rhoPhi_;
    const dimensionedScalar& rho1_;
    const dimensionedScalar& rho2_;

    fv::options& fvOptions_;

    IOdictionary thermalProperties_;

    dimensionedScalar cp1_;
    dimensionedScalar cp2_;
    dimensionedScalar kappa1_;
    dimensionedScalar kappa2_;
    dimensionedScalar Prt_;

    dimensionedScalar Tsat_;
    dimensionedScalar Lv_;
    dimensionedScalar Cevap_;
    dimensionedScalar Ccond_;

    volScalarField T_;
    volScalarField cp_;
    volScalarField kappa_;
    volScalarField rhoCp_;
    surfaceScalarField rhophicp_;

    volScalarField mDotEvap_;
    volScalarField mDotCond_;
    volScalarField mDot_;
    volScalarField latentSource_;
    volScalarField vDot_;

    volScalarField SuAlpha_;
    volScalarField SpAlpha_;

public:
    interIsoThermalEnergy
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const volScalarField& rho,
        const surfaceScalarField& rhoPhi,
        const dimensionedScalar& rho1,
        const dimensionedScalar& rho2,
        fv::options& fvOptions
    );

    void correct(); // Update properties
    void correctPhaseChange(); 
    void solve();

    const volScalarField& T() const; // Make T field public

    const volScalarField& mDot() const;
    const volScalarField& latentSource() const;
    const volScalarField& vDot() const;

    const volScalarField& SuAlpha() const;
    const volScalarField& SpAlpha() const;
};

}

#endif
