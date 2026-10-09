/*---------------------------------------------------------------------------*\
  Description
      Inline function file for LU decomposition of matrix.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/


inline double& FastChemistry::LUsolver::operator()
(
    const int i,
    const int j
)
{
    return v_[i*N + j];
}

inline double FastChemistry::LUsolver::operator() 
(
    const int i,
    const int j
)const
{
    return v_[i*N + j];
}

inline int& FastChemistry::LUsolver::pivotIndice
(
    const int i
)
{
    return this->pivotIndice_[i];
}

inline int FastChemistry::LUsolver::pivotIndice
(
    const int i
) const
{
    return this->pivotIndice_[i];
}