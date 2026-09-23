/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
     \\/     M anipulation  |
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


/*---------------------------------------------------------------------------*\
  Description
      Inline function for FastChemistryModel object

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class UnusedThermo>
inline label
Foam::FastChemistryModel<UnusedThermo>::nSpecie() const
{
    return nSpecie_;
}

template<class UnusedThermo>
inline Foam::PtrList<Foam::volScalarField>&
Foam::FastChemistryModel<UnusedThermo>::RR()
{
    return RR_;
}

template<class UnusedThermo>
inline Foam::label Foam::FastChemistryModel<UnusedThermo>::nReaction() const
{
    return GasReaction->n_Reactions;
}

template<class UnusedThermo>
inline std::unique_ptr<FastChemistry::OptReaction>& 
Foam::FastChemistryModel<UnusedThermo>::getReaction() const
{
    return GasReaction;
}


template<class UnusedThermo>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FastChemistryModel<UnusedThermo>::RR(const label i) const
{
    return RR_[i].internalField();
}


template<class UnusedThermo>
Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FastChemistryModel<UnusedThermo>::RR(const label i)
{
    return RR_[i].ref();
}


// ************************************************************************* //
