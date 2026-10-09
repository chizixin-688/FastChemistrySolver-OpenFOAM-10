/*---------------------------------------------------------------------------*\
  Description
      Header file for the block LU decomposition of a dense matrix.

      LUsolver decomposes an N x N matrix (stored row-major with a row stride
      of alignN = ceil(N/4)*4) into a unit lower-triangular part L and an
      upper-triangular part U using 4x4 diagonal blocks.  It is used by the
      fast-chemistry ODE solvers to solve the linear systems.

      The AVX-2 kernels are implemented in LUsolver.C.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <vector>

//---------------------------------
// 2. FastChemistry headers
//---------------------------------
#include "Macro.H"
#include "vec_math_avx2.H"

//---------------------------------
// 3. SIMD / AVX2 headers
//---------------------------------
#include <immintrin.h>

//=============================================================================//

#ifndef LUsolver_H
#define LUsolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FastChemistry
{
class LUsolver
{
    // Private data

        //- Matrix size.
        int N;

        //- Row stride of the storage (rounded up to a multiple of 4).
        int alignN;

        //- Data pointer (matrix is stored row-major, stride alignN).
        double* __restrict__ v_;

        //- Pivot index of each row (row permutation table).
        std::vector<int> pivotIndice_;
        
        //- Row start addresses; LUDecompose4 may swap entries of this array
        //  internally, so each block always sees its rows in the identity order.
        std::vector<double*> rowPtr;

        //- Reciprocal of the diagonal elements of U.
        std::vector<double> invD;

        //- Reciprocal of each row's maximum absolute value (row equilibration).
        std::vector<double> rowScales;

public:

    // Constructors

        //- Construct from externally owned data.
        LUsolver(double* externalData, int size);

        //- Disable default construction.
        LUsolver() = delete;

    // Destructor

        ~LUsolver();

    // Member functions

        // Access

            //- (i, j) element access (non-const).
            inline double& operator()(const int i,const int j);

            //- (i, j) element access (const).
            inline double operator() (const int i,const int j)const;

            //- pivotIndice access (non-const).
            inline int& pivotIndice(const int i);

            //- pivotIndice access (const).
            inline int pivotIndice(const int i) const;

        // Debug printing

            void printMatrix(double* A,int mRows, int nCol);

            void printMatrix();

            void printPivotIndice();

        // Block LU decomposition algorithm

            // using 4*4 block.
            void Block4LUDecompose();

        // Solve Ax=b.

            void xSolve_Serial(double* __restrict__ b);
            void xSolve(double* __restrict__ b);

        // Sub functions

            // LU decomposition of 4*4 block.
            void LUDecompose4(int k0);
            void LUDecompose4_2();
            void LUDecompose4_3();

            void forwardSubstitute4_0(int k0,int k1);
            void forwardSubstitute4_1(int k0,int k1);
            void forwardSubstitute4_2(int k0,int k1);
            void forwardSubstitute4_3(int k0,int k1);

            void backSubstitute4_0(int k0,int k1);
            void backSubstitute4_1(int k0,int k1);
            void backSubstitute4_2(int k0,int k1);
            void backSubstitute4_3(int k0,int k1);

            void permutation0(int k0,int k1);
            void permutation1(int k0,int k1);
            void permutation2(int k0,int k1);
            void permutation3(int k0,int k1);
            void permutation_2();
            void permutation_3();

            void UpdateL22U22_Vec_38(int k0,int k1);
            void scaleA0();
            void scaleA1();
            void scaleA2();
            void scaleA3();
            void scaleb(double* __restrict__ b);

        // Reassign a square matrix
        void ReAssign(double* externalData);

};


} // End namespace FastChemistry

#include "LUsolverI.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
