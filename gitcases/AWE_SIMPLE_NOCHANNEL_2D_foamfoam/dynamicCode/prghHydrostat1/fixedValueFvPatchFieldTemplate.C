/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 839b494fa341e145082a44e619d71ae1f6e4072b
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void prghHydrostat1_839b494fa341e145082a44e619d71ae1f6e4072b(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    prghHydrostat1FixedValueFvPatchScalarField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
prghHydrostat1FixedValueFvPatchScalarField::
prghHydrostat1FixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct prghHydrostat1 : patch/DimensionedField");
    }
}


Foam::
prghHydrostat1FixedValueFvPatchScalarField::
prghHydrostat1FixedValueFvPatchScalarField
(
    const prghHydrostat1FixedValueFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct prghHydrostat1 : patch/DimensionedField/mapper");
    }
}


Foam::
prghHydrostat1FixedValueFvPatchScalarField::
prghHydrostat1FixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct prghHydrostat1 : patch/dictionary");
    }
}


Foam::
prghHydrostat1FixedValueFvPatchScalarField::
prghHydrostat1FixedValueFvPatchScalarField
(
    const prghHydrostat1FixedValueFvPatchScalarField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct prghHydrostat1");
    }
}


Foam::
prghHydrostat1FixedValueFvPatchScalarField::
prghHydrostat1FixedValueFvPatchScalarField
(
    const prghHydrostat1FixedValueFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct prghHydrostat1 : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
prghHydrostat1FixedValueFvPatchScalarField::
~prghHydrostat1FixedValueFvPatchScalarField()
{
    if (false)
    {
        printMessage("Destroy prghHydrostat1");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
prghHydrostat1FixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs prghHydrostat1");
    }

//{{{ begin code
    #line 40 "/home/au779625/mySource/myPostdoc_AU/case/simple_2D/nochannel_FOAM_FOAM/AWE_SIMPLE_NOCHANNEL_2D_foamfoam/0/p_rgh.boundaryField.OUTLET"
const vectorField Cf(this->patch().Cf());
            const scalarField p_(patch().size(), Zero);
            scalarField prgh_(patch().size(), Zero);
            const scalarField rho_(this->patch().lookupPatchField<volScalarField, scalar>("rho"));
            const scalarField gh_(this->patch().lookupPatchField<volScalarField, scalar>("gh"));
            forAll(Cf,faceI)
            {

            prgh_[faceI] = p_[faceI]-rho_[faceI]*gh_[faceI];
            }
            operator==(prgh_);
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

