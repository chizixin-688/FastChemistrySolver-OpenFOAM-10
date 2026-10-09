/*---------------------------------------------------------------------------*\
  Description
      Header file for data block. This is used for MPI communication when
      dynamic load balance is used, for local Euler scheme

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

#ifndef dataBlock_H
#define dataBlock_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class dataBlock;

bool operator==(const dataBlock&, const dataBlock&);
bool operator!=(const dataBlock&, const dataBlock&);

Istream& operator >> (Istream& is, dataBlock& data);
Ostream& operator << (Ostream& os, const dataBlock& data);

/*---------------------------------------------------------------------------*\
                       Class dataBlock Declaration
\*---------------------------------------------------------------------------*/

class dataBlock
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

        //- the flow time step used in computation
        scalar deltaT;

        //- the chemical integration step in last computation
        scalar deltaTChem_;

public:

    // Constructors

        //- Default constructor
        dataBlock() = default;

        //- Construct from total number of species
        dataBlock(label& nSpecies_);

        //- Construct from Istream
        dataBlock(Istream& is);


    // Friend Operators

        inline friend bool operator==(const dataBlock&, const dataBlock&);
        inline friend bool operator!=(const dataBlock&, const dataBlock&);


    // IOstream operator

        friend Istream& operator>>(Istream&, dataBlock&);
        friend Ostream& operator<<(Ostream&, const dataBlock&);

};

//- inline operator
inline bool operator==(const dataBlock& data1, const dataBlock& data2)
{
    return
        data1.Y_ == data2.Y_ &&
        data1.T == data2.T &&
        data1.p == data2.p &&
        data1.CPUtime == data2.CPUtime &&
        data1.celli == data2.celli &&
        data1.deltaT == data2.deltaT &&
        data1.deltaTChem_ == data2.deltaTChem_;
}

inline bool operator!=(const dataBlock& data1, const dataBlock& data2)
{
    return !(data1==data2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
