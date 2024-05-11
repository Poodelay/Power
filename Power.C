/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "Power.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "cartesianCS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(Power, 0);
        addToRunTimeSelectionTable(functionObject, Power, dictionary);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::Power::createFiles()
{
    if (writeToFile() && !coeffFilePtr_)
    {
        coeffFilePtr_ = createFile("power");
        writeIntegratedHeader("Power", coeffFilePtr_());
    }
}

void Foam::functionObjects::Power::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, "Power coefficients");
    writeCommented(os, "Time");
    writeTabbed(os, "Pin");
    os  << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Power::Power
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool readFields
)
:
    forces(name, runTime, dict, false),
    power_({0,0,0}),
    coeffFilePtr_()
{
if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict, "liftDir", "dragDir");
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::Power::read(const dictionary& dict)
{
    forces::read(dict);
    
    dict.readEntry("magUInf", magUInf_);
    
    if (rhoName_ != "rhoInf")
    {
        dict.readEntry("rhoInf", rhoRef_);
    }
    
    dict.readEntry("lRef", lRef_);
    dict.readEntry("Aref", Aref_);
    
    return true;
}


bool Foam::functionObjects::Power::execute()
{
    forces::calcForcesMoment();

    createFiles();
    
    const volScalarField& p = lookupObject<volScalarField>(pName_);
    
    const volVectorField& U = lookupObject<volVectorField>(UName_);

    const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

    tmp<volSymmTensorField> tdevRhoReff = forces::devRhoReff();
    const volSymmTensorField::Boundary& devRhoReffb
        = tdevRhoReff().boundaryField();
        
    scalar pRef = pRef_/rho(p);

    power_ = {0,0,0};
    
    
    for (const label patchi : patchSet_)
    {
        vectorField Uc = U.boundaryField()[patchi]; //U.boundaryField()[patchi].patchInternalField();

        vectorField fT(Sfb[patchi] & devRhoReffb[patchi]);

        vectorField fN
        (
            rho(p)*Sfb[patchi]*(p.boundaryField()[patchi] - pRef)
        );
        
        vectorField f = fT + fN;
        
        power_[0] += gSum(Uc.component(0) * f.component(0));
        
        power_[1] += gSum(Uc.component(1) * f.component(1));
        
        power_[2] += gSum(Uc.component(2) * f.component(2));
    }
    
    //Pstream::listCombineGather(power_, plusEqOp<scalar>());
    //Pstream::listCombineScatter(power_);
    
    const scalar coeff = 1.0 / (0.5  * magUInf_ * magUInf_ * magUInf_ * Aref_);
    
    Log << "        coeff    : " << coeff << nl;     
    
    power_[0] *= coeff;
    power_[1] *= coeff;
    power_[2] *= coeff;
    
    if (writeToFile())
    {
        writeCurrentTime(coeffFilePtr_());
        coeffFilePtr_() << tab << power_.component(0) << tab << power_.component(1) << tab << power_.component(2) << endl;
    }
    
    Log << "        Pin      : " << power_ << nl;
    
    return true;
}


bool Foam::functionObjects::Power::write()
{
    return true;
}


// ************************************************************************* //
