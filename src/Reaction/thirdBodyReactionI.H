inline double FastChemistry::OptReaction::dMdC
(
    const unsigned int ri,
    const unsigned int cj
) const noexcept
{
    return this->ThirdBodyFactor1D[ri*this->AlignSpecies+cj];
}

