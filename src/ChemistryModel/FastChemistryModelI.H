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
    return reaction.n_Reactions;
}

template<class UnusedThermo>
inline FastChemistry::OptReaction& Foam::FastChemistryModel<UnusedThermo>::getReaction() const
{
    return reaction;
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


template<class UnusedThermo>
double Foam::FastChemistryModel<UnusedThermo>::hsum4(__m256d v) const
{
    __m128d lo = _mm256_castpd256_pd128(v);
    __m128d hi = _mm256_extractf128_pd(v,1);
    lo = _mm_add_pd(lo,hi);               // [v0+v2, v1+v3]
    __m128d sh = _mm_unpackhi_pd(lo,lo);  // [v1+v3, v1+v3]
    __m128d s = _mm_add_sd(lo, sh);       // [sum, …]
    return _mm_cvtsd_f64(s);
}

template<class UnusedThermo>
__m256d Foam::FastChemistryModel<UnusedThermo>::hsum4x4
(
    __m256d sum0,
    __m256d sum1,
    __m256d sum2,
    __m256d sum3
) const
{
    __m256d t0 = _mm256_hadd_pd(sum0, sum1);  // [a0+a1, b0+b1, a2+a3, b2+b3]
    __m256d t1 = _mm256_hadd_pd(sum2, sum3);  // [c0+c1, d0+d1, c2+c3, d2+d3]
    __m256d t2 = _mm256_permute4x64_pd (t0, _MM_SHUFFLE(3,1,2,0)); // [a0+a1, a2+a3, b0+b1,  b2+b3]
    __m256d t3 = _mm256_permute4x64_pd(t1, _MM_SHUFFLE(3,1,2,0)); // [c0+c1, c2+c3, d0+d1, d2+d3]
    __m256d t4 = _mm256_hadd_pd(t2, t3); // [a0+a1+a2+a3, c0+c1+c2+c3, b0+b1+b2+b3, d0+d1+d2+d3]
    __m256d t5 = _mm256_permute4x64_pd(t4, _MM_SHUFFLE(3,1,2,0));
    return t5;
}

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::transpose4x4_pd
(
    __m256d& v0,
    __m256d& v1,
    __m256d& v2,
    __m256d& v3
) const
{
    __m256d t0 = _mm256_unpacklo_pd(v0, v1); // a0 b0 a1 b1
    __m256d t1 = _mm256_unpackhi_pd(v0, v1); // a2 b2 a3 b3
    __m256d t2 = _mm256_unpacklo_pd(v2, v3); // c0 d0 c1 d1
    __m256d t3 = _mm256_unpackhi_pd(v2, v3); // c2 d2 c3 d3

    v0 = _mm256_permute2f128_pd(t0, t2, 0x20); // [a0 b0 c0 d0]
    v1 = _mm256_permute2f128_pd(t1, t3, 0x20); // [a1 b1 c1 d1]
    v2 = _mm256_permute2f128_pd(t0, t2, 0x31); // [a2 b2 c2 d2]
    v3 = _mm256_permute2f128_pd(t1, t3, 0x31); // [a3 b3 c3 d3]
}

template<class UnusedThermo>
double Foam::FastChemistryModel<UnusedThermo>::get_elem0(__m256d vec) const
{
    return _mm256_cvtsd_f64(vec); 
}
template<class UnusedThermo>
double Foam::FastChemistryModel<UnusedThermo>::get_elem1(__m256d vec) const
{
    __m128d low = _mm256_castpd256_pd128(vec); 
    return _mm_cvtsd_f64(_mm_unpackhi_pd(low, low));
}
template<class UnusedThermo>
double Foam::FastChemistryModel<UnusedThermo>::get_elem2(__m256d vec) const
{
    __m128d high = _mm256_extractf128_pd(vec, 1); 
    return _mm_cvtsd_f64(high);
}
template<class UnusedThermo>
double Foam::FastChemistryModel<UnusedThermo>::get_elem3(__m256d vec) const
{
    __m128d high = _mm256_extractf128_pd(vec, 1); 
    return _mm_cvtsd_f64(_mm_unpackhi_pd(high, high));
}

// ************************************************************************* //
