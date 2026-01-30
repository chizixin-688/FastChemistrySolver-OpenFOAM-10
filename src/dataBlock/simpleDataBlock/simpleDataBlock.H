/*---------------------------------------------------------------------------*\
  Description
      Header file for data block. This is used for MPI communication when
      dynamic load balance is used, for transient scheme

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/


//=============================================================================//

//---------------------------------
// 1. OpenFOAM library headers
//---------------------------------
#include "fvCFD.H"
#include "Ostream.H"
#include "Istream.H"

//=============================================================================//

#ifndef simpleDataBlock_H
#define simpleDataBlock_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class simpleDataBlock;

bool operator==(const simpleDataBlock&, const simpleDataBlock&);
bool operator!=(const simpleDataBlock&, const simpleDataBlock&);

Istream& operator >> (Istream& is, simpleDataBlock& data);
Ostream& operator << (Ostream& os, const simpleDataBlock& data);

/*---------------------------------------------------------------------------*\
                       Class simpleDataBlock Declaration
\*---------------------------------------------------------------------------*/

class simpleDataBlock
{

public:

    // Public Data
        
        //- a list of species mass fraction for chemistry reaction
        List<scalar> Y_;

        //- temperature [K]
        scalar T;

        //- absolute pressure [Pa]
        scalar p;

        //- CPU time for chemistry integration
        int64_t CPUtime;

        //- cell id
        label celli;

        //- the chemical integration step in last computation
        scalar deltaTChem_;

public:

    // Constructors

        //- Default constructor
        simpleDataBlock() = default;

        //- Construct from total number of species
        simpleDataBlock(label& nSpecies_);

        //- Construct from Istream
        simpleDataBlock(Istream& is);


    // Friend Operators

        inline friend bool operator==(const simpleDataBlock&, const simpleDataBlock&);
        inline friend bool operator!=(const simpleDataBlock&, const simpleDataBlock&);


    // IOstream operator

        friend Istream& operator>>(Istream&, simpleDataBlock&);
        friend Ostream& operator<<(Ostream&, const simpleDataBlock&);

};

//- inline operator
inline bool operator==(const simpleDataBlock& data1, const simpleDataBlock& data2)
{
    return
        data1.Y_ == data2.Y_ &&
        data1.T == data2.T &&
        data1.p == data2.p &&
        data1.CPUtime == data2.CPUtime &&
        data1.celli == data2.celli &&
        data1.deltaTChem_ == data2.deltaTChem_;
}

inline bool operator!=(const simpleDataBlock& data1, const simpleDataBlock& data2)
{
    return !(data1==data2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
