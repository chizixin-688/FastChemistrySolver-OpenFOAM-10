/*---------------------------------------------------------------------------*\
  Description
      Header file for data block. This is used for MPI communication when
      dynamic load balance is used, for transient scheme

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/
//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "simpleDataBlock.H"

//=============================================================================//



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleDataBlock::simpleDataBlock(label& nSpecies_):
Y_(nSpecies_),
T(300.0),
p(101325),
CPUtime(0),
celli(-1),
deltaTChem_(0)
{}

Foam::simpleDataBlock::simpleDataBlock(Istream& is)
{
    this->Y_.clear();
    is >> this->Y_;
    is >> this->T;
    is >> this->p;
    is >> this->CPUtime;
    is >> this->celli;
    is >> this->deltaTChem_;
}

// * * * * * * * * * * * * * * * IOstream Operator  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator >>(Istream& is, simpleDataBlock& data)
{
    is >> data.Y_;
    is >> data.T;
    is >> data.p;
    is >> data.CPUtime;
    is >> data.celli;
    is >> data.deltaTChem_;
    return is;
}

Foam::Ostream& Foam::operator <<(Ostream& os, const simpleDataBlock& data)
{
    os << data.Y_;
    os << data.T;
    os << data.p;
    os << data.CPUtime;
    os << data.celli;
    os << data.deltaTChem_;
    return os;
}






// ************************************************************************* //
