/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
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
    Foam::surfaceWriters::numpyWriter

Description
    A surfaceWriter for outputting in numpy format.
    The conceptual layout is similar to the boundaryData writer.

    \verbatim
    formatOptions
    {
        numpy
        {
            normal      false;
        }
    }
    \endverbatim

    Format options:
    \table
        Property    | Description                           | Required | Default
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
        normal      | Write face area normal in output      | no  | false
    \endtable

    Typical way of working:
    - use a sampledSurface of type 'patch' (to sample a patch):
    \verbatim
    surfaces
    {
        type            surfaces;
        fields          ( p );
        surfaceFormat   numpy;
        formatOptions
        {
            numpy
            {
                normal  yes;
                fieldLevel
                {
                    p       1e5;    // Absolute -> gauge [Pa]
                }
                fieldScale
                {
                    "p.*"   0.01;   // [Pa] -> [mbar]
                }
            }
        }
        surfaces
        {
            outlet
            {
                type            patch;
                patches         (outlet);
                interpolate     false;
            }
        }
    }
    \endverbatim

    \heading Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    where the geometry is written as:
    \verbatim
    rootdir
    `-- surfaceName
        `-- "points.npy"
    \endverbatim

    and field data:
    \verbatim
    rootdir
    `-- surfaceName
        |-- "points.npy"
        `-- timeName
            `-- field.npy
    \endverbatim

SourceFiles
    numpySurfaceWriter.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_numpyWriter_H
#define Foam_surfaceWriters_numpyWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                       Class numpyWriter Declaration
\*---------------------------------------------------------------------------*/

class numpyWriter
:
    public surfaceWriter
{
    // Private Data

        //- Output face area normal
        const bool writeNormal_;


    // Private Member Functions

        //- Write serial surface geometry to "points" file.
        void serialWriteGeometry
        (
            const fileName& directory,
            const meshedSurf& surf
        );

        //- Templated write field operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("numpy");


    // Constructors

        //- Default construct
        numpyWriter();

        //- Construct with some output options
        explicit numpyWriter(const dictionary& options);

        //- Construct from components
        numpyWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        numpyWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~numpyWriter() = default;


    // Member Functions

        //- Write surface geometry to file.
        virtual fileName write(); // override

        declareSurfaceWriterWriteMethod(label);
        declareSurfaceWriterWriteMethod(scalar);
        declareSurfaceWriterWriteMethod(vector);
        declareSurfaceWriterWriteMethod(sphericalTensor);
        declareSurfaceWriterWriteMethod(symmTensor);
        declareSurfaceWriterWriteMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
