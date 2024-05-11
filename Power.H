/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

Class
    Foam::functionObjects::UpdateBC

Group
    grpForcesFunctionObjects

Description
 
See also
    Foam::functionObject
    Foam::functionObjects::timeControl
    Foam::functionObjects::forces

SourceFiles
    UpdateBC.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_UpdateBC_H
#define functionObjects_UpdateBC_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class UpdateBC Declaration
\*---------------------------------------------------------------------------*/

class Power
:
    public forces
{
    // Private data
        
        // input Power
        
        vector power_;

    // File streams
    
        scalar magUInf_;


    // Reference scales

        //- Reference length [m]
        scalar lRef_;

        //- Reference area [m^2]
        scalar Aref_;

        //- Integrated coefficients
        autoPtr<OFstream> coeffFilePtr_;

    // Private Member Functions

        //- No copy construct
        Power(const Power&) = delete;

        //- No copy assignment
        void operator=(const Power&) = delete;


protected:

    // Protected Member Functions

        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeader(const word& header, Ostream& os) const;

public:

    //- Runtime type information
    TypeName("Power");


    // Constructors

        //- Construct from Time and dictionary
        Power
        (
            const word& name,
            const Time& runTime,
            const dictionary&,
            const bool readFields = true
        );


    //- Destructor
    virtual ~Power() = default;


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write the forces
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
