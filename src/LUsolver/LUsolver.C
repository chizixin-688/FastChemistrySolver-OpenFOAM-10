/*---------------------------------------------------------------------------*\
  Description
      Implementation of the block 4x4 LU decomposition kernels used by
      LUsolver.  The matrix is stored row-major with a row stride of
      alignN = ceil(N/4)*4.

      Main decomposition driver:  Block4LUDecompose()
        - scaleA0 .. scaleA3        : row equilibration (Remain-specific)
        - LUDecompose4              : in-place 4x4 block LU with pivoting
        - forwardSubstitute4_0..3   : forward elimination with the unit lower
                                      triangular block L (Remain-specific)
        - permutation0..3           : undo the block pivoting in the data
        - backSubstitute4_0..3      : solve with the upper triangular block U
                                      (Remain-specific)
        - UpdateL22U22_Vec_38       : trailing block update
      Solving Ax = b after the decomposition:
        - xSolve / xSolve_Serial 

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <iostream>
#include <cmath>
#include <iomanip>
#include <utility>

//---------------------------------
// 2. FastChemistry headers
//---------------------------------
#include "LUsolver.H"

//---------------------------------
// 3. SIMD / AVX2 headers
//---------------------------------
#include <immintrin.h>

//=============================================================================//

FastChemistry::LUsolver::LUsolver(double* externalData, int size)
{

    this->N = size;
    this->alignN = ((size+3)/4)*4;
    this->v_ = externalData;


    //this->pivotIndice_ = new int[N];
    this->pivotIndice_.resize(N);
    for(int i = 0; i< N; i++)
    {
        this->pivotIndice_[i] = i;
    }
    rowPtr.resize(N);
    for(int i = 0; i <N;i++)
    {
        rowPtr[i] = &v_[i*alignN];
    }
    invD.resize(N);
    this->rowScales.resize(N);
}

FastChemistry::LUsolver::~LUsolver()
{
    //delete[] this->pivotIndice_;
    //this->pivotIndice_ = nullptr;
    this->v_ = nullptr;
}


// Print an arbitrary array as a matrix (debug helper).
void FastChemistry::LUsolver::printMatrix(double*A,int mRows, int nCol)
{
    for(int i = 0;i<mRows;i++)
    {
        for(int j =0;j<nCol-1;j++)
        {
            std::cout<<std::setprecision(16)<<A[i*nCol+j]<<" ";
        }
        std::cout<<std::setprecision(16)<<A[i*nCol+nCol-1]<<std::endl;
    }
}

// Print the stored matrix (debug helper).
void FastChemistry::LUsolver::printMatrix()
{
    for(int i = 0;i<this->N;i++)
    {
        for(int j =0;j<this->N-1;j++)
        {
            std::cout<<this->v_[i*this->N+j]<<" ";
        }
        std::cout<<this->v_[i*this->N+this->N-1]<<std::endl;
    }
}

// Print the pivot indices (debug helper).
void FastChemistry::LUsolver::printPivotIndice()
{
    for(int i = 0;i<this->N;i++)
    {
        std::cout<<this->pivotIndice_[i]<<std::endl;
    }
}

// LU decomposition of the 4x4 diagonal block at column k0 (partial pivoting).
void FastChemistry::LUsolver::LUDecompose4
(
    int k0
)
{

// step 1
    int iMax = 0;
    double temp = 0;
    double* __restrict__ rowk0 = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1 = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2 = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3 = &v_[(k0+3)*alignN];


    {
        if(temp<std::fabs(rowk0[k0]))
        {
            iMax = 0;
            temp = std::fabs(rowk0[k0]);
        }
        if(temp<std::fabs(rowk1[k0]))
        {
            iMax = 1;
            temp = std::fabs(rowk1[k0]);
        }
        if(temp<std::fabs(rowk2[k0]))
        {
            iMax = 2;
            temp = std::fabs(rowk2[k0]);
        }
        if(temp<std::fabs(rowk3[k0]))
        {
            iMax = 3;
            temp = std::fabs(rowk3[k0]);
        }
    }

    if(iMax!=0)
    {
        std::swap(rowPtr[k0], rowPtr[k0+iMax]);
        this->pivotIndice_[0+k0] = iMax+k0;
        rowk0 = rowPtr[k0+0];
        rowk1 = rowPtr[k0+1];
        rowk2 = rowPtr[k0+2];
        rowk3 = rowPtr[k0+3];
    }
    if(rowk0[k0]==0)
    {
        rowk0[k0] = FastChemistry::LuLimiter;
    }
    double rU00 = 1.0/rowk0[k0];


    rowk1[k0] = rowk1[k0]*rU00;
    rowk2[k0] = rowk2[k0]*rU00;
    rowk3[k0] = rowk3[k0]*rU00;

// step 2

    rowk1[k0+1] = rowk1[k0+1] - rowk1[k0]*rowk0[k0+1];
    rowk2[k0+1] = rowk2[k0+1] - rowk2[k0]*rowk0[k0+1];
    rowk3[k0+1] = rowk3[k0+1] - rowk3[k0]*rowk0[k0+1];


    temp = std::fabs(rowk1[k0+1]);
    iMax = 1;

    if (temp<std::fabs(rowk2[k0+1]))
    {
        iMax = 2;
        temp = std::fabs(rowk2[k0+1]);
    }

    if (temp<std::fabs(rowk3[k0+1]))
    {
        iMax = 3;
        temp = std::fabs(rowk3[k0+1]);
    }

    if(iMax!=1)
    {
        std::swap(rowPtr[k0+1], rowPtr[k0+iMax]);
        rowk0 = rowPtr[k0+0];
        rowk1 = rowPtr[k0+1];
        rowk2 = rowPtr[k0+2];
        rowk3 = rowPtr[k0+3];
        this->pivotIndice_[1+k0] = iMax+k0;
    }

    if(rowk1[k0+1]==0)
    {
        rowk1[k0+1] = FastChemistry::LuLimiter;
    }

    rowk1[k0+2] = rowk1[k0+2] - rowk1[k0+0]*rowk0[k0+2];
    rowk1[k0+3] = rowk1[k0+3] - rowk1[k0+0]*rowk0[k0+3];

    double rU11 = 1.0/rowk1[k0+1];

    rowk2[k0+1] = rowk2[k0+1]*rU11;
    rowk3[k0+1] = rowk3[k0+1]*rU11;

// step 3

    rowk2[k0+2] = rowk2[k0+2] - rowk2[k0+0]*rowk0[k0+2] - rowk2[k0+1]*rowk1[k0+2];
    rowk3[k0+2] = rowk3[k0+2] - rowk3[k0+0]*rowk0[k0+2] - rowk3[k0+1]*rowk1[k0+2];

    if (std::fabs(rowk2[k0+2])<std::fabs(rowk3[k0+2]))
    {
        std::swap(rowPtr[k0+2], rowPtr[k0+3]);
        rowk0 = rowPtr[k0+0];
        rowk1 = rowPtr[k0+1];
        rowk2 = rowPtr[k0+2];
        rowk3 = rowPtr[k0+3];

        this->pivotIndice_[2+k0] ++;
    }
    if(rowk2[k0+2]==0)
    {
        rowk2[k0+2] = FastChemistry::LuLimiter;
    }

    rowk2[k0+3] =  rowk2[k0+3] - rowk2[k0+0]*rowk0[k0+3] - rowk2[k0+1]*rowk1[k0+3];

    double rU22 = 1.0/rowk2[k0+2];
    rowk3[k0+2] =  rowk3[k0+2] * rU22;


// step 4
    rowk3[k0+3] =  rowk3[k0+3] - rowk3[k0+0]*rowk0[k0+3] - rowk3[k0+1]*rowk1[k0+3] - rowk3[k0+2]*rowk2[k0+3];
    if(rowk3[k0+3]==0)
    {
        rowk3[k0+3] = FastChemistry::LuLimiter;
    }
    double rU33 = 1.0/rowk3[k0+3];
    this->invD[k0+0] = rU00;
    this->invD[k0+1] = rU11;
    this->invD[k0+2] = rU22;
    this->invD[k0+3] = rU33;
}

// LU decomposition of the trailing 2x2 block (rows/columns N-2, N-1).
void FastChemistry::LUsolver::LUDecompose4_2
(
)
{
    double* __restrict__ rowN2 = &v_[(N-2)*alignN];
    double* __restrict__ rowN1 = &v_[(N-1)*alignN];
    //[a b]
    //[c d]
    double a = rowN2[N-2];
    double b = rowN2[N-1];
    double c = rowN1[N-2];
    double d = rowN1[N-1];

    if( std::fabs(a) >= std::fabs(c))
    {
        if(a==0)
        {
            a = FastChemistry::LuLimiter;
        }
        double inva = 1.0/a;
        rowN1[N-2] = c*inva;
        rowN1[N-1] = d - c*b*inva;
        this->invD[N-2] = inva;
        if(rowN1[N-1]==0)
        {
            rowN1[N-1] = FastChemistry::LuLimiter;
        }
        this->invD[N-1] = 1.0/rowN1[N-1];
    }
    else if(std::fabs(a) < std::fabs(c))
    {
        if(c==0)
        {
            c = FastChemistry::LuLimiter;
        }
        double invc = 1.0/c;
        rowN2[N-2] = c;
        rowN2[N-1] = d;
        rowN1[N-2] = a*invc;
        rowN1[N-1] = b - a*d*invc;
        this->pivotIndice_[N-2] = N-1;
        this->invD[N-2] = invc;
        if(rowN1[N-1]==0)
        {
            rowN1[N-1] = FastChemistry::LuLimiter;
        }
        this->invD[N-1] = 1.0/rowN1[N-1];
    }
}

// LU decomposition of the trailing 3x3 block (rows/columns N-3 .. N-1).
void FastChemistry::LUsolver::LUDecompose4_3
(
)
{
    double Array0[4];
    {
        int k0 = this->N - 3;

// step 1
        int iMax = 0;
        double temp = 0;

        for(int i = 0; i < 3; i ++)
        {
            if(temp<std::fabs(rowPtr[k0+i][k0]))
            {
                iMax = i;
                temp = std::fabs(rowPtr[k0+i][k0]);
            }
        }
        double* __restrict__ rowk0 = &v_[(k0+0)*alignN];
        double* __restrict__ rowk1 = &v_[(k0+1)*alignN];
        double* __restrict__ rowk2 = &v_[(k0+2)*alignN];

        if(iMax!=0)
        {
            double* __restrict__ rowi = &v_[(k0+iMax)*alignN];
            Array0[0] = rowk0[k0+0];
            Array0[1] = rowk0[k0+1];
            Array0[2] = rowk0[k0+2];

            rowk0[k0+0] = rowi[k0+0];
            rowk0[k0+1] = rowi[k0+1];
            rowk0[k0+2] = rowi[k0+2];

            rowi[k0+0] = Array0[0];
            rowi[k0+1] = Array0[1];
            rowi[k0+2] = Array0[2];
            this->pivotIndice_[0+k0] = iMax+k0;
        }

        if(rowk0[k0+0]==0)
        {
            rowk0[k0+0] = FastChemistry::LuLimiter;
        }
        double rU00 = 1.0/rowk0[k0+0];

        rowk1[k0+0] = rowk1[k0+0]*rU00;
        rowk2[k0+0] = rowk2[k0+0]*rU00;

// step 2
        rowk1[k0+1] = rowk1[k0+1] - rowk1[k0+0]*rowk0[k0+1];
        rowk2[k0+1] = rowk2[k0+1] - rowk2[k0+0]*rowk0[k0+1];

        temp = std::fabs(rowk1[k0+1]);
        iMax = 1;

        if (temp<std::fabs(rowk2[k0+1]))
        {
            iMax = 2;
            temp = rowk2[k0+1];

            Array0[0] = rowk1[k0+0];
            Array0[1] = rowk1[k0+1];
            Array0[2] = rowk1[k0+2];

            rowk1[k0+0] = rowk2[k0+0];
            rowk1[k0+1] = rowk2[k0+1];
            rowk1[k0+2] = rowk2[k0+2];

            rowk2[k0+0] = Array0[0];
            rowk2[k0+1] = Array0[1];
            rowk2[k0+2] = Array0[2];

            this->pivotIndice_[1+k0] = 2+k0;
        }

        if(rowk1[k0+1]==0)
        {
            rowk1[k0+1] = FastChemistry::LuLimiter;
        }

        rowk1[k0+2] = rowk1[k0+2] - rowk1[k0+0]*rowk0[k0+2];
        double rU11 = 1.0/rowk1[k0+1];
        rowk2[k0+1] = rowk2[k0+1]*rU11;
        rowk2[k0+2] =  rowk2[k0+2] - rowk2[k0+0]*rowk0[k0+2] - rowk2[k0+1]*rowk1[k0+2];

        if(rowk2[k0+2]==0)
        {
            rowk2[k0+2] = FastChemistry::LuLimiter;
        }
        double rU22 = 1.0/rowk2[k0+2];
        this->invD[N-3] = rU00;
        this->invD[N-2] = rU11;
        this->invD[N-1] = rU22;
    }

}

// Forward substitution with the unit lower-triangular 4x4 block L:
//   rowk1 -= L10 * rowk0
//   rowk2 -= L20 * rowk0 + L21 * rowk1
//   rowk3 -= L30 * rowk0 + L31 * rowk1 + L32 * rowk2
// The four Remain variants only differ in how the (N-k1)%8 tail columns are
// handled.
void FastChemistry::LUsolver::forwardSubstitute4_0
(
    int k0,
    int k1
)
{
    double* __restrict__ rowk0 = rowPtr[k0+0];
    double* __restrict__ rowk1 = rowPtr[k0+1];
    double* __restrict__ rowk2 = rowPtr[k0+2];
    double* __restrict__ rowk3 = rowPtr[k0+3];

    const double L10 = rowk1[k0+0];
    const double L20 = rowk2[k0+0];
    const double L30 = rowk3[k0+0];
    const double L21 = rowk2[k0+1];
    const double L31 = rowk3[k0+1];
    const double L32 = rowk3[k0+2];
    __m256d L10v = _mm256_set1_pd(L10);
    __m256d L20v = _mm256_set1_pd(L20);
    __m256d L21v = _mm256_set1_pd(L21);
    __m256d L30v = _mm256_set1_pd(L30);
    __m256d L31v = _mm256_set1_pd(L31);
    __m256d L32v = _mm256_set1_pd(L32);


    int remain8 = (this->N-k1)%8;
    for(int i = k1; i < this->N-remain8; i=i+8)
    {
        __m256d U00v = load256d(&rowk0[i+0]);
        __m256d U04v = load256d(&rowk0[i+4]);

        __m256d U10v = load256d(&rowk1[i+0]);
        __m256d U14v = load256d(&rowk1[i+4]);
        U10v = fnmadd256d(L10v,U00v,U10v);
        U14v = fnmadd256d(L10v,U04v,U14v);
        store256d(&rowk1[i+0],U10v);
        store256d(&rowk1[i+4],U14v);

        __m256d U20v = load256d(&rowk2[i+0]);
        __m256d U24v = load256d(&rowk2[i+4]);
        U20v = fnmadd256d(L20v,U00v,U20v);
        U20v = fnmadd256d(L21v,U10v,U20v);
        U24v = fnmadd256d(L20v,U04v,U24v);
        U24v = fnmadd256d(L21v,U14v,U24v);
        store256d(&rowk2[i+0],U20v);
        store256d(&rowk2[i+4],U24v);

        __m256d U30v = load256d(&rowk3[i+0]);
        __m256d U34v = load256d(&rowk3[i+4]);
        U30v = fnmadd256d(L30v,U00v,U30v);
        U30v = fnmadd256d(L31v,U10v,U30v);
        U30v = fnmadd256d(L32v,U20v,U30v);
        U34v = fnmadd256d(L30v,U04v,U34v);
        U34v = fnmadd256d(L31v,U14v,U34v);
        U34v = fnmadd256d(L32v,U24v,U34v);
        store256d(&rowk3[i+0],U30v);
        store256d(&rowk3[i+4],U34v);
    }
    if(remain8==4)
    {
        int i = N-4;
        __m256d U00v = load256d(&rowk0[i+0]);
        __m256d U10v = load256d(&rowk1[i+0]);
        U10v = fnmadd256d(L10v,U00v,U10v);
        store256d(&rowk1[i+0],U10v);

        __m256d U20v = load256d(&rowk2[i+0]);
        U20v = fnmadd256d(L20v,U00v,U20v);
        U20v = fnmadd256d(L21v,U10v,U20v);
        store256d(&rowk2[i+0],U20v);

        __m256d U30v = load256d(&rowk3[i+0]);
        U30v = fnmadd256d(L30v,U00v,U30v);
        U30v = fnmadd256d(L31v,U10v,U30v);
        U30v = fnmadd256d(L32v,U20v,U30v);
        store256d(&rowk3[i+0],U30v);

    }
}

// Remain == 1 variant: tail is 4 columns at N-5 plus the scalar element N-1.
void FastChemistry::LUsolver::forwardSubstitute4_1
(
    int k0,
    int k1
)
{
    double* __restrict__ rowk0 = rowPtr[k0+0];
    double* __restrict__ rowk1 = rowPtr[k0+1];
    double* __restrict__ rowk2 = rowPtr[k0+2];
    double* __restrict__ rowk3 = rowPtr[k0+3];

    const double L10 = rowk1[k0+0];
    const double L20 = rowk2[k0+0];
    const double L30 = rowk3[k0+0];
    const double L21 = rowk2[k0+1];
    const double L31 = rowk3[k0+1];
    const double L32 = rowk3[k0+2];
    __m256d L10v = _mm256_set1_pd(L10);
    __m256d L20v = _mm256_set1_pd(L20);
    __m256d L21v = _mm256_set1_pd(L21);
    __m256d L30v = _mm256_set1_pd(L30);
    __m256d L31v = _mm256_set1_pd(L31);
    __m256d L32v = _mm256_set1_pd(L32);

    int remain8 = (this->N-k1)%8;
    for(int i = k1; i < this->N-remain8; i=i+8)
    {
        __m256d U00v = load256d(&rowk0[i+0]);
        __m256d U04v = load256d(&rowk0[i+4]);

        __m256d U10v = load256d(&rowk1[i+0]);
        __m256d U14v = load256d(&rowk1[i+4]);
        U10v = fnmadd256d(L10v,U00v,U10v);
        U14v = fnmadd256d(L10v,U04v,U14v);
        store256d(&rowk1[i+0],U10v);
        store256d(&rowk1[i+4],U14v);

        __m256d U20v = load256d(&rowk2[i+0]);
        __m256d U24v = load256d(&rowk2[i+4]);
        U20v = fnmadd256d(L20v,U00v,U20v);
        U20v = fnmadd256d(L21v,U10v,U20v);
        U24v = fnmadd256d(L20v,U04v,U24v);
        U24v = fnmadd256d(L21v,U14v,U24v);
        store256d(&rowk2[i+0],U20v);
        store256d(&rowk2[i+4],U24v);

        __m256d U30v = load256d(&rowk3[i+0]);
        __m256d U34v = load256d(&rowk3[i+4]);
        U30v = fnmadd256d(L30v,U00v,U30v);
        U30v = fnmadd256d(L31v,U10v,U30v);
        U30v = fnmadd256d(L32v,U20v,U30v);
        U34v = fnmadd256d(L30v,U04v,U34v);
        U34v = fnmadd256d(L31v,U14v,U34v);
        U34v = fnmadd256d(L32v,U24v,U34v);
        store256d(&rowk3[i+0],U30v);
        store256d(&rowk3[i+4],U34v);
    }
    if(remain8==5)
    {
        int i = N-5;
        __m256d U00v = load256d(&rowk0[i+0]);
        __m256d U10v = load256d(&rowk1[i+0]);
        U10v = fnmadd256d(L10v,U00v,U10v);
        store256d(&rowk1[i+0],U10v);

        __m256d U20v = load256d(&rowk2[i+0]);
        U20v = fnmadd256d(L20v,U00v,U20v);
        U20v = fnmadd256d(L21v,U10v,U20v);
        store256d(&rowk2[i+0],U20v);

        __m256d U30v = load256d(&rowk3[i+0]);
        U30v = fnmadd256d(L30v,U00v,U30v);
        U30v = fnmadd256d(L31v,U10v,U30v);
        U30v = fnmadd256d(L32v,U20v,U30v);
        store256d(&rowk3[i+0],U30v);

    }
    {
        {
            rowk1[N-1] = rowk1[N-1]
                         - L10*rowk0[N-1];
            rowk2[N-1] = rowk2[N-1]
                         - L20*rowk0[N-1]
                         - L21*rowk1[N-1];
            rowk3[N-1] = rowk3[N-1]
                         - L30*rowk0[N-1]
                         - L31*rowk1[N-1]
                         - L32*rowk2[N-1];
        }

    }
}

// Remain == 2 variant: tail is 4 columns at N-6 plus 2 columns at N-2.
void FastChemistry::LUsolver::forwardSubstitute4_2
(
    int k0,
    int k1
)
{
    double* __restrict__ rowk0 = rowPtr[k0+0];
    double* __restrict__ rowk1 = rowPtr[k0+1];
    double* __restrict__ rowk2 = rowPtr[k0+2];
    double* __restrict__ rowk3 = rowPtr[k0+3];

    const double L10 = rowk1[k0+0];
    const double L20 = rowk2[k0+0];
    const double L30 = rowk3[k0+0];
    const double L21 = rowk2[k0+1];
    const double L31 = rowk3[k0+1];
    const double L32 = rowk3[k0+2];
    __m256d L10v = _mm256_set1_pd(L10);
    __m256d L20v = _mm256_set1_pd(L20);
    __m256d L21v = _mm256_set1_pd(L21);
    __m256d L30v = _mm256_set1_pd(L30);
    __m256d L31v = _mm256_set1_pd(L31);
    __m256d L32v = _mm256_set1_pd(L32);
    int remain8 = (this->N-k1)%8;
    for(int i = k1; i < this->N-remain8; i=i+8)
    {
        __m256d U00v = load256d(&rowk0[i+0]);
        __m256d U04v = load256d(&rowk0[i+4]);

        __m256d U10v = load256d(&rowk1[i+0]);
        __m256d U14v = load256d(&rowk1[i+4]);
        U10v = fnmadd256d(L10v,U00v,U10v);
        U14v = fnmadd256d(L10v,U04v,U14v);
        store256d(&rowk1[i+0],U10v);
        store256d(&rowk1[i+4],U14v);

        __m256d U20v = load256d(&rowk2[i+0]);
        __m256d U24v = load256d(&rowk2[i+4]);
        U20v = fnmadd256d(L20v,U00v,U20v);
        U20v = fnmadd256d(L21v,U10v,U20v);
        U24v = fnmadd256d(L20v,U04v,U24v);
        U24v = fnmadd256d(L21v,U14v,U24v);
        store256d(&rowk2[i+0],U20v);
        store256d(&rowk2[i+4],U24v);

        __m256d U30v = load256d(&rowk3[i+0]);
        __m256d U34v = load256d(&rowk3[i+4]);
        U30v = fnmadd256d(L30v,U00v,U30v);
        U30v = fnmadd256d(L31v,U10v,U30v);
        U30v = fnmadd256d(L32v,U20v,U30v);
        U34v = fnmadd256d(L30v,U04v,U34v);
        U34v = fnmadd256d(L31v,U14v,U34v);
        U34v = fnmadd256d(L32v,U24v,U34v);
        store256d(&rowk3[i+0],U30v);
        store256d(&rowk3[i+4],U34v);
    }
    if(remain8==6)
    {
        int i = N-6;
        __m256d U00v = load256d(&rowk0[i+0]);
        __m256d U10v = load256d(&rowk1[i+0]);
        U10v = fnmadd256d(L10v,U00v,U10v);
        store256d(&rowk1[i+0],U10v);

        __m256d U20v = load256d(&rowk2[i+0]);
        U20v = fnmadd256d(L20v,U00v,U20v);
        U20v = fnmadd256d(L21v,U10v,U20v);
        store256d(&rowk2[i+0],U20v);

        __m256d U30v = load256d(&rowk3[i+0]);
        U30v = fnmadd256d(L30v,U00v,U30v);
        U30v = fnmadd256d(L31v,U10v,U30v);
        U30v = fnmadd256d(L32v,U20v,U30v);
        store256d(&rowk3[i+0],U30v);

    }
    {
        int i = N-2;
        __m128d U00v = load128d(&rowk0[i+0]);
        __m128d U10v = load128d(&rowk1[i+0]);
        U10v = fnmadd128d(get01(L10v),U00v,U10v);
        store128d(&rowk1[i+0],U10v);

        __m128d U20v = load128d(&rowk2[i+0]);
        U20v = fnmadd128d(get01(L20v),U00v,U20v);
        U20v = fnmadd128d(get01(L21v),U10v,U20v);
        store128d(&rowk2[i+0],U20v);

        __m128d U30v = load128d(&rowk3[i+0]);
        U30v = fnmadd128d(get01(L30v),U00v,U30v);
        U30v = fnmadd128d(get01(L31v),U10v,U30v);
        U30v = fnmadd128d(get01(L32v),U20v,U30v);
        store128d(&rowk3[i+0],U30v);
    }
}


// Remain == 3 variant: tail is 4 columns at N-7 plus 3 columns at N-3.
void FastChemistry::LUsolver::forwardSubstitute4_3
(
    int k0,
    int k1
)
{
    double* __restrict__ rowk0 = rowPtr[k0+0];
    double* __restrict__ rowk1 = rowPtr[k0+1];
    double* __restrict__ rowk2 = rowPtr[k0+2];
    double* __restrict__ rowk3 = rowPtr[k0+3];

    const double L10 = rowk1[k0+0];
    const double L20 = rowk2[k0+0];
    const double L30 = rowk3[k0+0];
    const double L21 = rowk2[k0+1];
    const double L31 = rowk3[k0+1];
    const double L32 = rowk3[k0+2];
    __m256d L10v = _mm256_set1_pd(L10);
    __m256d L20v = _mm256_set1_pd(L20);
    __m256d L21v = _mm256_set1_pd(L21);
    __m256d L30v = _mm256_set1_pd(L30);
    __m256d L31v = _mm256_set1_pd(L31);
    __m256d L32v = _mm256_set1_pd(L32);


    int remain8 = (this->N-k1)%8;
    for(int i = k1; i < this->N-remain8; i=i+8)
    {
        __m256d U00v = load256d(&rowk0[i+0]);
        __m256d U04v = load256d(&rowk0[i+4]);

        __m256d U10v = load256d(&rowk1[i+0]);
        __m256d U14v = load256d(&rowk1[i+4]);
        U10v = fnmadd256d(L10v,U00v,U10v);
        U14v = fnmadd256d(L10v,U04v,U14v);
        store256d(&rowk1[i+0],U10v);
        store256d(&rowk1[i+4],U14v);

        __m256d U20v = load256d(&rowk2[i+0]);
        __m256d U24v = load256d(&rowk2[i+4]);
        U20v = fnmadd256d(L20v,U00v,U20v);
        U20v = fnmadd256d(L21v,U10v,U20v);
        U24v = fnmadd256d(L20v,U04v,U24v);
        U24v = fnmadd256d(L21v,U14v,U24v);
        store256d(&rowk2[i+0],U20v);
        store256d(&rowk2[i+4],U24v);

        __m256d U30v = load256d(&rowk3[i+0]);
        __m256d U34v = load256d(&rowk3[i+4]);
        U30v = fnmadd256d(L30v,U00v,U30v);
        U30v = fnmadd256d(L31v,U10v,U30v);
        U30v = fnmadd256d(L32v,U20v,U30v);
        U34v = fnmadd256d(L30v,U04v,U34v);
        U34v = fnmadd256d(L31v,U14v,U34v);
        U34v = fnmadd256d(L32v,U24v,U34v);
        store256d(&rowk3[i+0],U30v);
        store256d(&rowk3[i+4],U34v);
    }
    if(remain8==7)
    {
        int i = N-7;
        __m256d U00v = load256d(&rowk0[i+0]);
        __m256d U10v = load256d(&rowk1[i+0]);
        U10v = fnmadd256d(L10v,U00v,U10v);
        store256d(&rowk1[i+0],U10v);

        __m256d U20v = load256d(&rowk2[i+0]);
        U20v = fnmadd256d(L20v,U00v,U20v);
        U20v = fnmadd256d(L21v,U10v,U20v);
        store256d(&rowk2[i+0],U20v);

        __m256d U30v = load256d(&rowk3[i+0]);
        U30v = fnmadd256d(L30v,U00v,U30v);
        U30v = fnmadd256d(L31v,U10v,U30v);
        U30v = fnmadd256d(L32v,U20v,U30v);
        store256d(&rowk3[i+0],U30v);

    }
    {
        int i = N-3;

        __m256d U00v = load256d012(&rowk0[i+0]);
        __m256d U10v = load256d012(&rowk1[i+0]);
        U10v = fnmadd256d(L10v,U00v,U10v);
        store256d012(&rowk1[i+0],U10v);

        __m256d U20v = load256d012(&rowk2[i+0]);
        U20v = fnmadd256d(L20v,U00v,U20v);
        U20v = fnmadd256d(L21v,U10v,U20v);
        store256d012(&rowk2[i+0],U20v);

        __m256d U30v = load256d012(&rowk3[i+0]);
        U30v = fnmadd256d(L30v,U00v,U30v);
        U30v = fnmadd256d(L31v,U10v,U30v);
        U30v = fnmadd256d(L32v,U20v,U30v);
        store256d012(&rowk3[i+0],U30v);
    }
}



// Back substitution: for every group of rows below the diagonal block, solve
// the 4x4 upper-triangular system with the block U (columns k0..k0+3):
//   x0 = l0 / U00
//   x1 = (l1 - U01*x0) / U11
//   x2 = (l2 - U02*x0 - U12*x1) / U22
//   x3 = (l3 - U03*x0 - U13*x1 - U23*x2) / U33
// and write the solution back into columns k0..k0+3.
void FastChemistry::LUsolver::backSubstitute4_0
(
    int k0,
    int k1
)
{
    double* __restrict__ rowk0 = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1 = rowk0 + alignN;
    double* __restrict__ rowk2 = rowk1 + alignN;

    const double invU00 = this->invD[k0+0];
    const double invU11 = this->invD[k0+1];
    const double invU22 = this->invD[k0+2];
    const double invU33 = this->invD[k0+3];
    const double U01 = rowk0[k0+1];
    const double U02 = rowk0[k0+2];
    const double U03 = rowk0[k0+3];
    const double U12 = rowk1[k0+2];
    const double U13 = rowk1[k0+3];
    const double U23 = rowk2[k0+3];
    __m256d invU00v = _mm256_set1_pd(invU00);
    __m256d invU11v = _mm256_set1_pd(invU11);
    __m256d invU22v = _mm256_set1_pd(invU22);
    __m256d invU33v = _mm256_set1_pd(invU33);
    __m256d U01v = _mm256_set1_pd(U01);
    __m256d U02v = _mm256_set1_pd(U02);
    __m256d U03v = _mm256_set1_pd(U03);
    __m256d U12v = _mm256_set1_pd(U12);
    __m256d U13v = _mm256_set1_pd(U13);
    __m256d U23v = _mm256_set1_pd(U23);


    for(int i = k1; i<this->N; i = i+4)
    {
        double* __restrict__ rowi0 = &v_[(i+0)*alignN];
        double* __restrict__ rowi1 = &v_[(i+1)*alignN];
        double* __restrict__ rowi2 = &v_[(i+2)*alignN];
        double* __restrict__ rowi3 = &v_[(i+3)*alignN];

        __m256d L00v = load256d(&rowi0[k0+0]);
        __m256d L01v = load256d(&rowi1[k0+0]);
        __m256d L02v = load256d(&rowi2[k0+0]);
        __m256d L03v = load256d(&rowi3[k0+0]);
        transpose4x4_pd(L00v,L01v,L02v,L03v);

        L00v = _mm256_mul_pd(L00v,invU00v);

        L01v = fnmadd256d(L00v,U01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = fnmadd256d(L00v,U02v,L02v);
        L02v = fnmadd256d(L01v,U12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = fnmadd256d(L00v,U03v,L03v);
        L03v = fnmadd256d(L01v,U13v,L03v);
        L03v = fnmadd256d(L02v,U23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);
    }
}


// Remain == 1 variant: tail is 4 rows at N-5 plus the scalar row N-1.
void FastChemistry::LUsolver::backSubstitute4_1
(
    int k0,
    int k1
)
{
    double* __restrict__ rowk0 = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1 = rowk0 + alignN;
    double* __restrict__ rowk2 = rowk1 + alignN;

    const double invU00 = this->invD[k0+0];
    const double invU11 = this->invD[k0+1];
    const double invU22 = this->invD[k0+2];
    const double invU33 = this->invD[k0+3];
    const double U01 = rowk0[k0+1];
    const double U02 = rowk0[k0+2];
    const double U03 = rowk0[k0+3];
    const double U12 = rowk1[k0+2];
    const double U13 = rowk1[k0+3];
    const double U23 = rowk2[k0+3];
    __m256d invU00v = _mm256_set1_pd(invU00);
    __m256d invU11v = _mm256_set1_pd(invU11);
    __m256d invU22v = _mm256_set1_pd(invU22);
    __m256d invU33v = _mm256_set1_pd(invU33);
    __m256d U01v = _mm256_set1_pd(U01);
    __m256d U02v = _mm256_set1_pd(U02);
    __m256d U03v = _mm256_set1_pd(U03);
    __m256d U12v = _mm256_set1_pd(U12);
    __m256d U13v = _mm256_set1_pd(U13);
    __m256d U23v = _mm256_set1_pd(U23);

    int remain8 = (N-k1)%8;
    for(int i = k1; i < this->N-remain8; i=i+8)
    {
        double* __restrict__ rowi0 = &v_[(i+0)*alignN];
        double* __restrict__ rowi1 = &v_[(i+1)*alignN];
        double* __restrict__ rowi2 = &v_[(i+2)*alignN];
        double* __restrict__ rowi3 = &v_[(i+3)*alignN];
        __m256d L00v = load256d(&rowi0[k0+0]);
        __m256d L01v = load256d(&rowi1[k0+0]);
        __m256d L02v = load256d(&rowi2[k0+0]);
        __m256d L03v = load256d(&rowi3[k0+0]);
        transpose4x4_pd(L00v,L01v,L02v,L03v);

        L00v = _mm256_mul_pd(L00v,invU00v);

        L01v = fnmadd256d(L00v,U01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = fnmadd256d(L00v,U02v,L02v);
        L02v = fnmadd256d(L01v,U12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = fnmadd256d(L00v,U03v,L03v);
        L03v = fnmadd256d(L01v,U13v,L03v);
        L03v = fnmadd256d(L02v,U23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);
        double* __restrict__ rowi4 = &v_[(i+4)*alignN];
        double* __restrict__ rowi5 = &v_[(i+5)*alignN];
        double* __restrict__ rowi6 = &v_[(i+6)*alignN];
        double* __restrict__ rowi7 = &v_[(i+7)*alignN];
        __m256d L40v = load256d(&rowi4[k0+0]);
        __m256d L41v = load256d(&rowi5[k0+0]);
        __m256d L42v = load256d(&rowi6[k0+0]);
        __m256d L43v = load256d(&rowi7[k0+0]);
        transpose4x4_pd(L40v,L41v,L42v,L43v);

        L40v = _mm256_mul_pd(L40v,invU00v);

        L41v = fnmadd256d(L40v,U01v,L41v);
        L41v = _mm256_mul_pd(L41v,invU11v);

        L42v = fnmadd256d(L40v,U02v,L42v);
        L42v = fnmadd256d(L41v,U12v,L42v);
        L42v = _mm256_mul_pd(L42v,invU22v);

        L43v = fnmadd256d(L40v,U03v,L43v);
        L43v = fnmadd256d(L41v,U13v,L43v);
        L43v = fnmadd256d(L42v,U23v,L43v);
        L43v = _mm256_mul_pd(L43v,invU33v);

        transpose4x4_pd(L40v,L41v,L42v,L43v);
        store256d(&rowi4[k0+0],L40v);
        store256d(&rowi5[k0+0],L41v);
        store256d(&rowi6[k0+0],L42v);
        store256d(&rowi7[k0+0],L43v);
    }
    if(remain8==5)
    {
        int i = N-5;
        double* __restrict__ rowi0 = &v_[(i+0)*alignN];
        double* __restrict__ rowi1 = &v_[(i+1)*alignN];
        double* __restrict__ rowi2 = &v_[(i+2)*alignN];
        double* __restrict__ rowi3 = &v_[(i+3)*alignN];

        __m256d L00v = load256d(&rowi0[k0+0]);
        __m256d L01v = load256d(&rowi1[k0+0]);
        __m256d L02v = load256d(&rowi2[k0+0]);
        __m256d L03v = load256d(&rowi3[k0+0]);
        transpose4x4_pd(L00v,L01v,L02v,L03v);

        L00v = _mm256_mul_pd(L00v,invU00v);

        L01v = fnmadd256d(L00v,U01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = fnmadd256d(L00v,U02v,L02v);
        L02v = fnmadd256d(L01v,U12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = fnmadd256d(L00v,U03v,L03v);
        L03v = fnmadd256d(L01v,U13v,L03v);
        L03v = fnmadd256d(L02v,U23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);
    }
    {
        double* __restrict__ rowN1 = &v_[(N-1)*alignN];
        double& La0 = rowN1[k0+0];
        La0 = La0*invU00;

        double& La1 = rowN1[k0+1];
        La1 = (La1 - La0*U01)*invU11;

        double& La2 = rowN1[k0+2];
        La2 = (La2 - La0*U02 - La1*U12)*invU22;

        double& La3 = rowN1[k0+3];
        La3 = (La3 - La0*U03 - La1*U13 - La2*U23)*invU33;
    }
}

// Remain == 2 variant: tail is 4 rows at N-6 plus 2 rows at N-2.
void FastChemistry::LUsolver::backSubstitute4_2
(
    int k0,
    int k1
)
{
    double* __restrict__ rowk0 = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1 = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2 = &v_[(k0+2)*alignN];

    const double invU00 = this->invD[k0+0];
    const double invU11 = this->invD[k0+1];
    const double invU22 = this->invD[k0+2];
    const double invU33 = this->invD[k0+3];
    const double U01 = rowk0[k0+1];
    const double U02 = rowk0[k0+2];
    const double U03 = rowk0[k0+3];
    const double U12 = rowk1[k0+2];
    const double U13 = rowk1[k0+3];
    const double U23 = rowk2[k0+3];
    __m256d invU00v = _mm256_set1_pd(invU00);
    __m256d invU11v = _mm256_set1_pd(invU11);
    __m256d invU22v = _mm256_set1_pd(invU22);
    __m256d invU33v = _mm256_set1_pd(invU33);
    __m256d U01v = _mm256_set1_pd(U01);
    __m256d U02v = _mm256_set1_pd(U02);
    __m256d U03v = _mm256_set1_pd(U03);
    __m256d U12v = _mm256_set1_pd(U12);
    __m256d U13v = _mm256_set1_pd(U13);
    __m256d U23v = _mm256_set1_pd(U23);
    int remain8 = (N-k1)%8;
    for(int i=k1; i<this->N-remain8; i=i+8)
    {
        double* __restrict__ rowi0 = &v_[(i+0)*alignN];
        double* __restrict__ rowi1 = &v_[(i+1)*alignN];
        double* __restrict__ rowi2 = &v_[(i+2)*alignN];
        double* __restrict__ rowi3 = &v_[(i+3)*alignN];
        __m256d L00v = load256d(&rowi0[k0+0]);
        __m256d L01v = load256d(&rowi1[k0+0]);
        __m256d L02v = load256d(&rowi2[k0+0]);
        __m256d L03v = load256d(&rowi3[k0+0]);

        transpose4x4_pd(L00v,L01v,L02v,L03v);

        L00v = _mm256_mul_pd(L00v,invU00v);

        L01v = fnmadd256d(L00v,U01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = fnmadd256d(L00v,U02v,L02v);
        L02v = fnmadd256d(L01v,U12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = fnmadd256d(L00v,U03v,L03v);
        L03v = fnmadd256d(L01v,U13v,L03v);
        L03v = fnmadd256d(L02v,U23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);

        double* __restrict__ rowi4 = &v_[(i+4)*alignN];
        double* __restrict__ rowi5 = &v_[(i+5)*alignN];
        double* __restrict__ rowi6 = &v_[(i+6)*alignN];
        double* __restrict__ rowi7 = &v_[(i+7)*alignN];
        __m256d L40v = load256d(&rowi4[k0+0]);
        __m256d L41v = load256d(&rowi5[k0+0]);
        __m256d L42v = load256d(&rowi6[k0+0]);
        __m256d L43v = load256d(&rowi7[k0+0]);

        transpose4x4_pd(L40v,L41v,L42v,L43v);

        L40v = _mm256_mul_pd(L40v,invU00v);

        L41v = fnmadd256d(L40v,U01v,L41v);
        L41v = _mm256_mul_pd(L41v,invU11v);

        L42v = fnmadd256d(L40v,U02v,L42v);
        L42v = fnmadd256d(L41v,U12v,L42v);
        L42v = _mm256_mul_pd(L42v,invU22v);

        L43v = fnmadd256d(L40v,U03v,L43v);
        L43v = fnmadd256d(L41v,U13v,L43v);
        L43v = fnmadd256d(L42v,U23v,L43v);
        L43v = _mm256_mul_pd(L43v,invU33v);

        transpose4x4_pd(L40v,L41v,L42v,L43v);
        store256d(&rowi4[k0+0],L40v);
        store256d(&rowi5[k0+0],L41v);
        store256d(&rowi6[k0+0],L42v);
        store256d(&rowi7[k0+0],L43v);
    }
    if(remain8==6)
    {
        int i = this->N-6;
        double* __restrict__ rowi0 = &v_[(i+0)*alignN];
        double* __restrict__ rowi1 = &v_[(i+1)*alignN];
        double* __restrict__ rowi2 = &v_[(i+2)*alignN];
        double* __restrict__ rowi3 = &v_[(i+3)*alignN];
        __m256d L00v = load256d(&rowi0[k0+0]);
        __m256d L01v = load256d(&rowi1[k0+0]);
        __m256d L02v = load256d(&rowi2[k0+0]);
        __m256d L03v = load256d(&rowi3[k0+0]);

        transpose4x4_pd(L00v,L01v,L02v,L03v);

        L00v = _mm256_mul_pd(L00v,invU00v);

        L01v = fnmadd256d(L00v,U01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = fnmadd256d(L00v,U02v,L02v);
        L02v = fnmadd256d(L01v,U12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = fnmadd256d(L00v,U03v,L03v);
        L03v = fnmadd256d(L01v,U13v,L03v);
        L03v = fnmadd256d(L02v,U23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);

    }
    {
        int i = this->N-2;
        double* __restrict__ rowi0 = &v_[(i+0)*alignN];
        double* __restrict__ rowi1 = &v_[(i+1)*alignN];

        __m256d L00v = load256d(&rowi0[k0+0]);//L00 L01 L02 L03
        __m256d L10v = load256d(&rowi1[k0+0]);//L10 L11 L12 L13

        __m128d v0 = _mm_setzero_pd();
        __m128d v1 = _mm_setzero_pd();
        __m128d v2 = _mm_setzero_pd();
        __m128d v3 = _mm_setzero_pd();
        transpose2x4_pd(L00v,L10v,v0,v1,v2,v3);

        v0 = _mm_mul_pd(v0,get01(invU00v));

        v1 = fnmadd128d(v0,get01(U01v),v1);
        v1 = _mm_mul_pd(v1,get01(invU11v));

        v2 = fnmadd128d(v0,get01(U02v),v2);
        v2 = fnmadd128d(v1,get01(U12v),v2);
        v2 = _mm_mul_pd(v2,get01(invU22v));

        v3 = fnmadd128d(v0,get01(U03v),v3);
        v3 = fnmadd128d(v1,get01(U13v),v3);
        v3 = fnmadd128d(v2,get01(U23v),v3);
        v3 = _mm_mul_pd(v3,get01(invU33v));

        transpose4x2_pd(L00v,L10v,v0,v1,v2,v3);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L10v);

    }
}
// Remain == 3 variant: tail is 4 rows at N-7 plus 3 rows at N-3.
void FastChemistry::LUsolver::backSubstitute4_3
(
    int k0,
    int k1
)
{
    double* __restrict__ rowk0 = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1 = rowk0 + alignN;
    double* __restrict__ rowk2 = rowk1 + alignN;

    const double invU00 = this->invD[k0+0];
    const double invU11 = this->invD[k0+1];
    const double invU22 = this->invD[k0+2];
    const double invU33 = this->invD[k0+3];
    const double U01 = rowk0[k0+1];
    const double U02 = rowk0[k0+2];
    const double U03 = rowk0[k0+3];
    const double U12 = rowk1[k0+2];
    const double U13 = rowk1[k0+3];
    const double U23 = rowk2[k0+3];
    __m256d invU00v = _mm256_set1_pd(invU00);
    __m256d invU11v = _mm256_set1_pd(invU11);
    __m256d invU22v = _mm256_set1_pd(invU22);
    __m256d invU33v = _mm256_set1_pd(invU33);
    __m256d U01v = _mm256_set1_pd(U01);
    __m256d U02v = _mm256_set1_pd(U02);
    __m256d U03v = _mm256_set1_pd(U03);
    __m256d U12v = _mm256_set1_pd(U12);
    __m256d U13v = _mm256_set1_pd(U13);
    __m256d U23v = _mm256_set1_pd(U23);
    int remain8 = (N-k1)%8;
    for(int i = k1; i < this->N-remain8; i=i+8)
    {
        double* __restrict__ rowi0 = &v_[(i+0)*alignN];
        double* __restrict__ rowi1 = rowi0+alignN;
        double* __restrict__ rowi2 = rowi1+alignN;
        double* __restrict__ rowi3 = rowi2+alignN;

        __m256d L00v = load256d(&rowi0[k0+0]);
        __m256d L01v = load256d(&rowi1[k0+0]);
        __m256d L02v = load256d(&rowi2[k0+0]);
        __m256d L03v = load256d(&rowi3[k0+0]);
        transpose4x4_pd(L00v,L01v,L02v,L03v);

        L00v = _mm256_mul_pd(L00v,invU00v);

        L01v = fnmadd256d(L00v,U01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = fnmadd256d(L00v,U02v,L02v);
        L02v = fnmadd256d(L01v,U12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = fnmadd256d(L00v,U03v,L03v);
        L03v = fnmadd256d(L01v,U13v,L03v);
        L03v = fnmadd256d(L02v,U23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);

        double* __restrict__ rowi4 = &v_[(i+4)*alignN];
        double* __restrict__ rowi5 = &v_[(i+5)*alignN];
        double* __restrict__ rowi6 = &v_[(i+6)*alignN];
        double* __restrict__ rowi7 = &v_[(i+7)*alignN];
        __m256d L40v = load256d(&rowi4[k0+0]);
        __m256d L41v = load256d(&rowi5[k0+0]);
        __m256d L42v = load256d(&rowi6[k0+0]);
        __m256d L43v = load256d(&rowi7[k0+0]);

        transpose4x4_pd(L40v,L41v,L42v,L43v);

        L40v = _mm256_mul_pd(L40v,invU00v);

        L41v = fnmadd256d(L40v,U01v,L41v);
        L41v = _mm256_mul_pd(L41v,invU11v);

        L42v = fnmadd256d(L40v,U02v,L42v);
        L42v = fnmadd256d(L41v,U12v,L42v);
        L42v = _mm256_mul_pd(L42v,invU22v);

        L43v = fnmadd256d(L40v,U03v,L43v);
        L43v = fnmadd256d(L41v,U13v,L43v);
        L43v = fnmadd256d(L42v,U23v,L43v);
        L43v = _mm256_mul_pd(L43v,invU33v);

        transpose4x4_pd(L40v,L41v,L42v,L43v);
        store256d(&rowi4[k0+0],L40v);
        store256d(&rowi5[k0+0],L41v);
        store256d(&rowi6[k0+0],L42v);
        store256d(&rowi7[k0+0],L43v);
    }
    if(remain8==7)
    {
        int i = N-7;
        double* __restrict__ rowi0 = &v_[(i+0)*alignN];
        double* __restrict__ rowi1 = &v_[(i+1)*alignN];
        double* __restrict__ rowi2 = &v_[(i+2)*alignN];
        double* __restrict__ rowi3 = &v_[(i+3)*alignN];

        __m256d L00v = load256d(&rowi0[k0+0]);
        __m256d L01v = load256d(&rowi1[k0+0]);
        __m256d L02v = load256d(&rowi2[k0+0]);
        __m256d L03v = load256d(&rowi3[k0+0]);
        transpose4x4_pd(L00v,L01v,L02v,L03v);

        L00v = _mm256_mul_pd(L00v,invU00v);

        L01v = fnmadd256d(L00v,U01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = fnmadd256d(L00v,U02v,L02v);
        L02v = fnmadd256d(L01v,U12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = fnmadd256d(L00v,U03v,L03v);
        L03v = fnmadd256d(L01v,U13v,L03v);
        L03v = fnmadd256d(L02v,U23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);
    }
    {
        int i = N-3;

        double* __restrict__ rowi0 = &v_[(i+0)*alignN];
        double* __restrict__ rowi1 = &v_[(i+1)*alignN];
        double* __restrict__ rowi2 = &v_[(i+2)*alignN];

        __m256d L00v = load256d(&rowi0[k0+0]);
        __m256d L01v = load256d(&rowi1[k0+0]);
        __m256d L02v = load256d(&rowi2[k0+0]);
        __m256d L03v = _mm256_setzero_pd();
        transpose4x4_pd(L00v,L01v,L02v,L03v);

        L00v = _mm256_mul_pd(L00v,invU00v);

        L01v = fnmadd256d(L00v,U01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = fnmadd256d(L00v,U02v,L02v);
        L02v = fnmadd256d(L01v,U12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = fnmadd256d(L00v,U03v,L03v);
        L03v = fnmadd256d(L01v,U13v,L03v);
        L03v = fnmadd256d(L02v,U23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
    }

}



// Trailing block update: A[k1:N][k1:N] -= L[k1:N][k0:k0+3] * U[k0:k0+3][k1:N].
void FastChemistry::LUsolver::UpdateL22U22_Vec_38
(
    int k0,
    int k1
)
{

    int rowRemain = (this->N-k1)%3;
    int colRemain = (this->alignN-k1)%8;
    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];
    for(int i = k1; i<this->N-rowRemain; i=i+3)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];

        __m256d L003v = load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);

        __m256d L203v = load256d(&rowi2[k0+0]);
        __m256d L20 = _mm256_permute4x64_pd(L203v, 0x00);
        __m256d L21 = _mm256_permute4x64_pd(L203v, 0x55);
        __m256d L22 = _mm256_permute4x64_pd(L203v, 0xAA);
        __m256d L23 = _mm256_permute4x64_pd(L203v, 0xFF);
        for(int j=k1; j<this->alignN-colRemain; j=j+8)
        {
            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14
            __m256d C = load256d(&rowi2[j]);        //15

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fnmadd_pd(L00,U003v,A);       //16
            B = _mm256_fnmadd_pd(L10,U003v,B);       //16
            C = _mm256_fnmadd_pd(L20,U003v,C);       //15

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fnmadd_pd(L01,U103v,A);       //16
            B = _mm256_fnmadd_pd(L11,U103v,B);       //16
            C = _mm256_fnmadd_pd(L21,U103v,C);       //15

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fnmadd_pd(L02,U203v,A);       //16
            B = _mm256_fnmadd_pd(L12,U203v,B);       //16
            C = _mm256_fnmadd_pd(L22,U203v,C);       //15

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fnmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fnmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14
            C = _mm256_fnmadd_pd(L23,U303v,C);       //13
            store256d(&rowi2[j], C);                //12


            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d E = load256d(&rowi1[j+4]);      //14
            __m256d F = load256d(&rowi2[j+4]);      //15

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fnmadd_pd(L00,U047v,D);       //11
            E = _mm256_fnmadd_pd(L10,U047v,E);       //10
            F = _mm256_fnmadd_pd(L20,U047v,F);

            __m256d U147v = load256d(&rowk1[j+4]);
            D = _mm256_fnmadd_pd(L01,U147v,D);
            E = _mm256_fnmadd_pd(L11,U147v,E);
            F = _mm256_fnmadd_pd(L21,U147v,F);

            __m256d U247v = load256d(&rowk2[j+4]);
            D = _mm256_fnmadd_pd(L02,U247v,D);
            E = _mm256_fnmadd_pd(L12,U247v,E);
            F = _mm256_fnmadd_pd(L22,U247v,F);

            __m256d U347v = load256d(&rowk3[j+4]);
            D = _mm256_fnmadd_pd(L03,U347v,D);
            store256d(&rowi0[j+4], D);
            E = _mm256_fnmadd_pd(L13,U347v,E);
            store256d(&rowi1[j+4], E);
            F = _mm256_fnmadd_pd(L23,U347v,F);
            store256d(&rowi2[j+4], F);
        }
        if(colRemain==4)
        {
            int j = this->alignN-colRemain;
            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14
            __m256d C = load256d(&rowi2[j]);        //15

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fnmadd_pd(L00,U003v,A);       //16
            B = _mm256_fnmadd_pd(L10,U003v,B);       //16
            C = _mm256_fnmadd_pd(L20,U003v,C);       //15

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fnmadd_pd(L01,U103v,A);       //16
            B = _mm256_fnmadd_pd(L11,U103v,B);       //16
            C = _mm256_fnmadd_pd(L21,U103v,C);       //15

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fnmadd_pd(L02,U203v,A);       //16
            B = _mm256_fnmadd_pd(L12,U203v,B);       //16
            C = _mm256_fnmadd_pd(L22,U203v,C);       //15

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fnmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fnmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14
            C = _mm256_fnmadd_pd(L23,U303v,C);       //13
            store256d(&rowi2[j], C);                //12
        }
    }
    if(rowRemain==2)
    {
        int i = this->N-2;
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];

        __m256d L003v = load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);

        for(int j=k1; j<this->alignN-colRemain; j=j+8)
        {
            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fnmadd_pd(L00,U003v,A);       //16
            B = _mm256_fnmadd_pd(L10,U003v,B);       //16

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fnmadd_pd(L01,U103v,A);       //16
            B = _mm256_fnmadd_pd(L11,U103v,B);       //16

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fnmadd_pd(L02,U203v,A);       //16
            B = _mm256_fnmadd_pd(L12,U203v,B);       //16

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fnmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fnmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14

            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d E = load256d(&rowi1[j+4]);      //14

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fnmadd_pd(L00,U047v,D);       //11
            E = _mm256_fnmadd_pd(L10,U047v,E);       //10

            __m256d U147v = load256d(&rowk1[j+4]);
            D = _mm256_fnmadd_pd(L01,U147v,D);
            E = _mm256_fnmadd_pd(L11,U147v,E);

            __m256d U247v = load256d(&rowk2[j+4]);
            D = _mm256_fnmadd_pd(L02,U247v,D);
            E = _mm256_fnmadd_pd(L12,U247v,E);

            __m256d U347v = load256d(&rowk3[j+4]);
            D = _mm256_fnmadd_pd(L03,U347v,D);
            store256d(&rowi0[j+4], D);
            E = _mm256_fnmadd_pd(L13,U347v,E);
            store256d(&rowi1[j+4], E);
        }
        if(colRemain==4)
        {
            int j = this->alignN-colRemain;
            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fnmadd_pd(L00,U003v,A);       //16
            B = _mm256_fnmadd_pd(L10,U003v,B);       //16

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fnmadd_pd(L01,U103v,A);       //16
            B = _mm256_fnmadd_pd(L11,U103v,B);       //16

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fnmadd_pd(L02,U203v,A);       //16
            B = _mm256_fnmadd_pd(L12,U203v,B);       //16

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fnmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fnmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14
        }
    }
    else if(rowRemain==1)
    {
        int i = this->N-1;
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];

        __m256d L003v = load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        for(int j=k1; j<this->alignN-colRemain; j=j+8)
        {
            __m256d A = load256d(&rowi0[j]);        //13
            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fnmadd_pd(L00,U003v,A);       //16

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fnmadd_pd(L01,U103v,A);       //16

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fnmadd_pd(L02,U203v,A);       //16

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fnmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15

            __m256d D = load256d(&rowi0[j+4]);      //13

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fnmadd_pd(L00,U047v,D);       //11

            __m256d U147v = load256d(&rowk1[j+4]);
            D = _mm256_fnmadd_pd(L01,U147v,D);

            __m256d U247v = load256d(&rowk2[j+4]);
            D = _mm256_fnmadd_pd(L02,U247v,D);

            __m256d U347v = load256d(&rowk3[j+4]);
            D = _mm256_fnmadd_pd(L03,U347v,D);
            store256d(&rowi0[j+4], D);
        }
        if(colRemain==4)
        {
            int j = this->alignN-colRemain;
            __m256d A = load256d(&rowi0[j]);        //13

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fnmadd_pd(L00,U003v,A);       //16

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fnmadd_pd(L01,U103v,A);       //16

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fnmadd_pd(L02,U203v,A);       //16

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fnmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
        }
    }
}



// Undo the block pivoting of block k0 by swapping the matrix rows back.
void FastChemistry::LUsolver::permutation0
(
    int k0,
    int k1
)
{
    (void)k1;  // k1 is not used by the permutation kernels.
    int remain16 = this->N%16;
    int remain8 = this->N%8;
    int remain4 = this->N%4;
    for(int j = 0; j<4;j++)
    {
        if(this->pivotIndice_[k0+j]==k0+j)
        {
            continue;
        }
        int iTarget = this->pivotIndice_[j+k0];
        double* __restrict__ rowiTarget = &v_[iTarget*alignN];
        double* __restrict__ rowj = &v_[(k0+j)*alignN];
        for(int i = 0; i<this->N-remain16; i=i+16)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d A1 = load256d(&rowiTarget[i+4]);
            __m256d A2 = load256d(&rowiTarget[i+8]);
            __m256d A3 = load256d(&rowiTarget[i+12]);
            __m256d B0 = load256d(&rowj[i+0]);
            __m256d B1 = load256d(&rowj[i+4]);
            __m256d B2 = load256d(&rowj[i+8]);
            __m256d B3 = load256d(&rowj[i+12]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowiTarget[i+4],B1);
            store256d(&rowiTarget[i+8],B2);
            store256d(&rowiTarget[i+12],B3);
            store256d(&rowj[i+0],A0);
            store256d(&rowj[i+4],A1);
            store256d(&rowj[i+8],A2);
            store256d(&rowj[i+12],A3);
        }
        for(int i = this->N-remain16; i<this->N-remain8; i=i+8)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d A1 = load256d(&rowiTarget[i+4]);
            __m256d B0 = load256d(&rowj[i+0]);
            __m256d B1 = load256d(&rowj[i+4]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowiTarget[i+4],B1);
            store256d(&rowj[i+0],A0);
            store256d(&rowj[i+4],A1);
        }
        for(int i = this->N-remain8; i<this->N-remain4; i=i+4)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d B0 = load256d(&rowj[i+0]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowj[i+0],A0);
        }
    }
}

// Remain == 1 variant of the pivot undo.
void FastChemistry::LUsolver::permutation1
(
    int k0,
    int k1
)
{
    (void)k1;  // k1 is not used by the permutation kernels.
    int remain16 = this->N%16;
    int remain8 = this->N%8;
    int remain4 = this->N%4;
    for(int j = 0; j<4;j++)
    {
        if(this->pivotIndice_[k0+j]==k0+j)
        {
            continue;
        }
        int iTarget = this->pivotIndice_[j+k0];
        double* __restrict__ rowiTarget = &v_[iTarget*alignN];
        double* __restrict__ rowj = &v_[(k0+j)*alignN];
        for(int i = 0; i<this->N-remain16; i=i+16)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d A1 = load256d(&rowiTarget[i+4]);
            __m256d A2 = load256d(&rowiTarget[i+8]);
            __m256d A3 = load256d(&rowiTarget[i+12]);
            __m256d B0 = load256d(&rowj[i+0]);
            __m256d B1 = load256d(&rowj[i+4]);
            __m256d B2 = load256d(&rowj[i+8]);
            __m256d B3 = load256d(&rowj[i+12]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowiTarget[i+4],B1);
            store256d(&rowiTarget[i+8],B2);
            store256d(&rowiTarget[i+12],B3);
            store256d(&rowj[i+0],A0);
            store256d(&rowj[i+4],A1);
            store256d(&rowj[i+8],A2);
            store256d(&rowj[i+12],A3);
        }
        for(int i = this->N-remain16; i<this->N-remain8; i=i+8)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d A1 = load256d(&rowiTarget[i+4]);
            __m256d B0 = load256d(&rowj[i+0]);
            __m256d B1 = load256d(&rowj[i+4]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowiTarget[i+4],B1);
            store256d(&rowj[i+0],A0);
            store256d(&rowj[i+4],A1);
        }
        for(int i = this->N-remain8; i<this->N-remain4; i=i+4)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d B0 = load256d(&rowj[i+0]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowj[i+0],A0);
        }
        /*for(int i = 0; i <= this->N-4; i=i+4)
        {
            __m256d A = load256d(&rowiTarget[i]);
            __m256d B = load256d(&rowj[i]);
            store256d(&rowiTarget[i],B);
            store256d(&rowj[i],A);
        }*/
        {
            int i = this->N-1;
            double A = rowiTarget[i];
            rowiTarget[i] = rowj[i];
            rowj[i] = A;
        }
    }
}

// Remain == 2 variant of the pivot undo.
void FastChemistry::LUsolver::permutation2
(
    int k0,
    int k1
)
{
    (void)k1;  // k1 is not used by the permutation kernels.
    int remain16 = this->N%16;
    int remain8 = this->N%8;
    int remain4 = this->N%4;
    for(int j = 0; j<4;j++)
    {
        if(this->pivotIndice_[k0+j]==k0+j)
        {
            continue;
        }
        int iTarget = this->pivotIndice_[j+k0];
        double* __restrict__ rowiTarget = &v_[iTarget*alignN];
        double* __restrict__ rowj = &v_[(k0+j)*alignN];
        for(int i = 0; i<this->N-remain16; i=i+16)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d A1 = load256d(&rowiTarget[i+4]);
            __m256d A2 = load256d(&rowiTarget[i+8]);
            __m256d A3 = load256d(&rowiTarget[i+12]);
            __m256d B0 = load256d(&rowj[i+0]);
            __m256d B1 = load256d(&rowj[i+4]);
            __m256d B2 = load256d(&rowj[i+8]);
            __m256d B3 = load256d(&rowj[i+12]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowiTarget[i+4],B1);
            store256d(&rowiTarget[i+8],B2);
            store256d(&rowiTarget[i+12],B3);
            store256d(&rowj[i+0],A0);
            store256d(&rowj[i+4],A1);
            store256d(&rowj[i+8],A2);
            store256d(&rowj[i+12],A3);
        }
        for(int i = this->N-remain16; i<this->N-remain8; i=i+8)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d A1 = load256d(&rowiTarget[i+4]);
            __m256d B0 = load256d(&rowj[i+0]);
            __m256d B1 = load256d(&rowj[i+4]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowiTarget[i+4],B1);
            store256d(&rowj[i+0],A0);
            store256d(&rowj[i+4],A1);
        }
        for(int i = this->N-remain8; i<this->N-remain4; i=i+4)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d B0 = load256d(&rowj[i+0]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowj[i+0],A0);
        }
        {
            int i = this->N-2;

            double A0 = rowiTarget[i+0];
            double A1 = rowiTarget[i+1];
            rowiTarget[i+0] = rowj[i+0];
            rowiTarget[i+1] = rowj[i+1];
            rowj[i+0] = A0;
            rowj[i+1] = A1;
        }
    }
}

// Remain == 3 variant of the pivot undo.
void FastChemistry::LUsolver::permutation3
(
    int k0,
    int k1
)
{
    (void)k1;  // k1 is not used by the permutation kernels.
    int remain16 = this->N%16;
    int remain8 = this->N%8;
    int remain4 = this->N%4;
    for(int j = 0; j<4;j++)
    {
        if(this->pivotIndice_[k0+j]==k0+j)
        {
            continue;
        }
        int iTarget = this->pivotIndice_[j+k0];
        double* __restrict__ rowiTarget = &v_[iTarget*alignN];
        double* __restrict__ rowj = &v_[(k0+j)*alignN];
        for(int i = 0; i<this->N-remain16; i=i+16)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d A1 = load256d(&rowiTarget[i+4]);
            __m256d A2 = load256d(&rowiTarget[i+8]);
            __m256d A3 = load256d(&rowiTarget[i+12]);
            __m256d B0 = load256d(&rowj[i+0]);
            __m256d B1 = load256d(&rowj[i+4]);
            __m256d B2 = load256d(&rowj[i+8]);
            __m256d B3 = load256d(&rowj[i+12]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowiTarget[i+4],B1);
            store256d(&rowiTarget[i+8],B2);
            store256d(&rowiTarget[i+12],B3);
            store256d(&rowj[i+0],A0);
            store256d(&rowj[i+4],A1);
            store256d(&rowj[i+8],A2);
            store256d(&rowj[i+12],A3);
        }
        for(int i = this->N-remain16; i<this->N-remain8; i=i+8)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d A1 = load256d(&rowiTarget[i+4]);
            __m256d B0 = load256d(&rowj[i+0]);
            __m256d B1 = load256d(&rowj[i+4]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowiTarget[i+4],B1);
            store256d(&rowj[i+0],A0);
            store256d(&rowj[i+4],A1);
        }
        for(int i = this->N-remain8; i<this->N-remain4; i=i+4)
        {
            __m256d A0 = load256d(&rowiTarget[i+0]);
            __m256d B0 = load256d(&rowj[i+0]);
            store256d(&rowiTarget[i+0],B0);
            store256d(&rowj[i+0],A0);
        }
        {
            int i = this->N-3;

            double A0 = rowiTarget[i+0];
            double A1 = rowiTarget[i+1];
            double A2 = rowiTarget[i+2];

            rowiTarget[i+0] = rowj[i+0];
            rowiTarget[i+1] = rowj[i+1];
            rowiTarget[i+2] = rowj[i+2];

            rowj[i+0] = A0;
            rowj[i+1] = A1;
            rowj[i+2] = A2;
        }
    }
}

// Undo the pivoting of the trailing 2 rows (N-2, N-1).
void FastChemistry::LUsolver::permutation_2
(
)
{
    double Array0[4];
    {
        if(this->pivotIndice_[N-2]==N-1)
        {
            for(int i = 0; i < this->N-2; i=i+4)
            {

                Array0[0] = v_[(N-1)*alignN+i+0];
                Array0[1] = v_[(N-1)*alignN+i+1];
                Array0[2] = v_[(N-1)*alignN+i+2];
                Array0[3] = v_[(N-1)*alignN+i+3];

                v_[(N-1)*alignN+i+0] = v_[(N-2)*alignN+i+0];
                v_[(N-1)*alignN+i+1] = v_[(N-2)*alignN+i+1];
                v_[(N-1)*alignN+i+2] = v_[(N-2)*alignN+i+2];
                v_[(N-1)*alignN+i+3] = v_[(N-2)*alignN+i+3];

                v_[(N-2)*alignN+i+0] = Array0[0];
                v_[(N-2)*alignN+i+1] = Array0[1];
                v_[(N-2)*alignN+i+2] = Array0[2];
                v_[(N-2)*alignN+i+3] = Array0[3];
            }
        }
    }
}

// Undo the pivoting of the trailing 3 rows (N-3, N-2, N-1).
void FastChemistry::LUsolver::permutation_3
(
)
{
    double Array0[4];
    {
        int k0 = this->N-3;
        for(int j = 0; j<3;j++)
        {
            if(this->pivotIndice_[k0+j]==k0+j)
            {
                continue;
            }

            for(int i = 0; i < k0; i=i+4)
            {
                int iTarget = this->pivotIndice_[j+k0];

                Array0[0] = v_[iTarget*alignN+(i+0)];
                Array0[1] = v_[iTarget*alignN+(i+1)];
                Array0[2] = v_[iTarget*alignN+(i+2)];
                Array0[3] = v_[iTarget*alignN+(i+3)];

                v_[iTarget*alignN+(i+0)] = v_[(k0+j)*alignN+i+0];
                v_[iTarget*alignN+(i+1)] = v_[(k0+j)*alignN+i+1];
                v_[iTarget*alignN+(i+2)] = v_[(k0+j)*alignN+i+2];
                v_[iTarget*alignN+(i+3)] = v_[(k0+j)*alignN+i+3];

                v_[(k0+j)*alignN+i+0] = Array0[0];
                v_[(k0+j)*alignN+i+1] = Array0[1];
                v_[(k0+j)*alignN+i+2] = Array0[2];
                v_[(k0+j)*alignN+i+3] = Array0[3];
            }
        }
    }

}


// Solve LUx = b (vectorized AVX-2 kernels).
void FastChemistry::LUsolver::xSolve
(
    double* __restrict__ b
)
{
    for(int i = 0; i < this->N; i++)
    {
        b[i] *= this->rowScales[i];
    }
    for(int j = 0; j<this->N;j++)
    {
        if(this->pivotIndice_[j]==j)
        {
            continue;
        }

        int jTarget = this->pivotIndice_[j];
        double temp = b[jTarget];
        b[jTarget] = b[j];
        b[j] = temp;
    }


    int remain = this->N%12;
    for(int i=0;i<this->N-remain;i=i+12)
    {

        const double* __restrict__ Li0 = &v_[(i+0)*alignN+0];
        const double* __restrict__ Li1 = &v_[(i+1)*alignN+0];
        const double* __restrict__ Li2 = &v_[(i+2)*alignN+0];
        const double* __restrict__ Li3 = &v_[(i+3)*alignN+0];
        const double* __restrict__ Li4 = &v_[(i+4)*alignN+0];
        const double* __restrict__ Li5 = &v_[(i+5)*alignN+0];
        const double* __restrict__ Li6 = &v_[(i+6)*alignN+0];
        const double* __restrict__ Li7 = &v_[(i+7)*alignN+0];
        const double* __restrict__ Li8 = &v_[(i+8)*alignN+0];
        const double* __restrict__ Li9 = &v_[(i+9)*alignN+0];
        const double* __restrict__ Li10 = &v_[(i+10)*alignN+0];
        const double* __restrict__ Li11 = &v_[(i+11)*alignN+0];



        __m256d yi0 = _mm256_setzero_pd();
        __m256d yi1 = _mm256_setzero_pd();
        __m256d yi2 = _mm256_setzero_pd();
        __m256d yi3 = _mm256_setzero_pd();
        __m256d yi4 = _mm256_setzero_pd();
        __m256d yi5 = _mm256_setzero_pd();
        __m256d yi6 = _mm256_setzero_pd();
        __m256d yi7 = _mm256_setzero_pd();
        __m256d yi8 = _mm256_setzero_pd();
        __m256d yi9 = _mm256_setzero_pd();
        __m256d yi10 = _mm256_setzero_pd();
        __m256d yi11 = _mm256_setzero_pd();
        for(int j = 0;j<i;j=j+12)
        {
            __m256d yjv = load256d(&b[j]);
            yi0 = fmadd256d(load256d(&Li0[j]),yjv,yi0);
            yi1 = fmadd256d(load256d(&Li1[j]),yjv,yi1);
            yi2 = fmadd256d(load256d(&Li2[j]),yjv,yi2);
            yi3 = fmadd256d(load256d(&Li3[j]),yjv,yi3);
            yi4 = fmadd256d(load256d(&Li4[j]),yjv,yi4);
            yi5 = fmadd256d(load256d(&Li5[j]),yjv,yi5);
            yi6 = fmadd256d(load256d(&Li6[j]),yjv,yi6);
            yi7 = fmadd256d(load256d(&Li7[j]),yjv,yi7);
            yi8 = fmadd256d(load256d(&Li8[j]),yjv,yi8);
            yi9 = fmadd256d(load256d(&Li9[j]),yjv,yi9);
            yi10 = fmadd256d(load256d(&Li10[j]),yjv,yi10);
            yi11 = fmadd256d(load256d(&Li11[j]),yjv,yi11);

            yjv = load256d(&b[j+4]);
            yi0 = fmadd256d(load256d(&Li0[j+4]),yjv,yi0);
            yi1 = fmadd256d(load256d(&Li1[j+4]),yjv,yi1);
            yi2 = fmadd256d(load256d(&Li2[j+4]),yjv,yi2);
            yi3 = fmadd256d(load256d(&Li3[j+4]),yjv,yi3);
            yi4 = fmadd256d(load256d(&Li4[j+4]),yjv,yi4);
            yi5 = fmadd256d(load256d(&Li5[j+4]),yjv,yi5);
            yi6 = fmadd256d(load256d(&Li6[j+4]),yjv,yi6);
            yi7 = fmadd256d(load256d(&Li7[j+4]),yjv,yi7);
            yi8 = fmadd256d(load256d(&Li8[j+4]),yjv,yi8);
            yi9 = fmadd256d(load256d(&Li9[j+4]),yjv,yi9);
            yi10 = fmadd256d(load256d(&Li10[j+4]),yjv,yi10);
            yi11 = fmadd256d(load256d(&Li11[j+4]),yjv,yi11);


            yjv = load256d(&b[j+8]);
            yi0 = fmadd256d(load256d(&Li0[j+8]),yjv,yi0);
            yi1 = fmadd256d(load256d(&Li1[j+8]),yjv,yi1);
            yi2 = fmadd256d(load256d(&Li2[j+8]),yjv,yi2);
            yi3 = fmadd256d(load256d(&Li3[j+8]),yjv,yi3);
            yi4 = fmadd256d(load256d(&Li4[j+8]),yjv,yi4);
            yi5 = fmadd256d(load256d(&Li5[j+8]),yjv,yi5);
            yi6 = fmadd256d(load256d(&Li6[j+8]),yjv,yi6);
            yi7 = fmadd256d(load256d(&Li7[j+8]),yjv,yi7);
            yi8 = fmadd256d(load256d(&Li8[j+8]),yjv,yi8);
            yi9 = fmadd256d(load256d(&Li9[j+8]),yjv,yi9);
            yi10 = fmadd256d(load256d(&Li10[j+8]),yjv,yi10);
            yi11 = fmadd256d(load256d(&Li11[j+8]),yjv,yi11);
        }
        //方案1
        
        __m256d zero = _mm256_setzero_pd();

        __m256d bi0v = load256d(&b[i+0]);

        __m256d bi0 = _mm256_blend_pd(bi0v, zero, 0b1110);  
        yi0 = sub256d(bi0,yi0);                                 

        __m256d bi1 = _mm256_blend_pd(bi0v, zero, 0b1101);  
        __m256d Li1_0 = _mm256_set1_pd(Li1[i+0]);
        yi1 = sub256d(bi1,yi1);
        yi1 = fnmadd256d(Li1_0, yi0, yi1);

        __m256d bi2 = _mm256_blend_pd(bi0v, zero, 0b1011);  
        yi2 = sub256d(bi2,yi2);
        __m256d Li2_0 = _mm256_set1_pd(Li2[i+0]);
        yi2 = fnmadd256d(Li2_0, yi0, yi2);
        __m256d Li2_1 = _mm256_set1_pd(Li2[i+1]);
        yi2 = fnmadd256d(Li2_1, yi1, yi2); 

        __m256d bi3 = _mm256_blend_pd(bi0v, zero, 0b0111);
        yi3 = sub256d(bi3,yi3);
        __m256d Li3_0 = _mm256_set1_pd(Li3[i+0]);
        yi3 = fnmadd256d(Li3_0, yi0, yi3);
        __m256d Li3_1 = _mm256_set1_pd(Li3[i+1]);
        yi3 = fnmadd256d(Li3_1, yi1, yi3);
        __m256d Li3_2 = _mm256_set1_pd(Li3[i+2]);
        yi3 = fnmadd256d(Li3_2, yi2, yi3);

        __m256d yi0v = hsum4x4(yi0,yi1,yi2,yi3);
        store256d(&b[i],yi0v);

        // update 2nd 4x4 block and 3rd 4x4 block
        yi4 = fmadd256d(load256d(&Li4[i]),yi0v,yi4);
        yi5 = fmadd256d(load256d(&Li5[i]),yi0v,yi5);
        yi6 = fmadd256d(load256d(&Li6[i]),yi0v,yi6);
        yi7 = fmadd256d(load256d(&Li7[i]),yi0v,yi7);
        yi8 = fmadd256d(load256d(&Li8[i]),yi0v,yi8);
        yi9 = fmadd256d(load256d(&Li9[i]),yi0v,yi9);
        yi10 = fmadd256d(load256d(&Li10[i]),yi0v,yi10);
        yi11 = fmadd256d(load256d(&Li11[i]),yi0v,yi11);

        // Update 2nd 4x4 diagonal block
        __m256d bi4v = load256d(&b[i+4]);                   

        __m256d bi4 = _mm256_blend_pd(bi4v, zero, 0b1110);  
        yi4 = sub256d(bi4,yi4);                                  

        __m256d bi5 = _mm256_blend_pd(bi4v, zero, 0b1101);  
        __m256d Li5_4 = _mm256_set1_pd(Li5[i+4]);
        yi5 = sub256d(bi5,yi5);
        yi5 = fnmadd256d(Li5_4, yi4, yi5);

        __m256d bi6 = _mm256_blend_pd(bi4v, zero, 0b1011);
        __m256d Li6_4 = _mm256_set1_pd(Li6[i+4]);
        __m256d Li6_5 = _mm256_set1_pd(Li6[i+5]);
        yi6 = sub256d(bi6,yi6);
        yi6 = fnmadd256d(Li6_4, yi4, yi6);
        yi6 = fnmadd256d(Li6_5, yi5, yi6);

        __m256d bi7 = _mm256_blend_pd(bi4v, zero, 0b0111);
        __m256d Li7_4 = _mm256_set1_pd(Li7[i+4]);
        __m256d Li7_5 = _mm256_set1_pd(Li7[i+5]);
        __m256d Li7_6 = _mm256_set1_pd(Li7[i+6]);
        yi7 = sub256d(bi7,yi7);
        yi7 = fnmadd256d(Li7_4, yi4, yi7);
        yi7 = fnmadd256d(Li7_5, yi5, yi7);
        yi7 = fnmadd256d(Li7_6, yi6, yi7);

        __m256d yi4v = hsum4x4(yi4,yi5,yi6,yi7);
        store256d(&b[i+4],yi4v);

        yi8 = fmadd256d(load256d(&Li8[i+4]),yi4v,yi8);
        yi9 = fmadd256d(load256d(&Li9[i+4]),yi4v,yi9);
        yi10 = fmadd256d(load256d(&Li10[i+4]),yi4v,yi10);
        yi11 = fmadd256d(load256d(&Li11[i+4]),yi4v,yi11);

        __m256d bi8v = load256d(&b[i+8]);

        __m256d bi8 = _mm256_blend_pd(bi8v, zero, 0b1110);
        yi8 = sub256d(bi8,yi8);

        __m256d bi9 = _mm256_blend_pd(bi8v, zero, 0b1101);
        __m256d Li9_8 = _mm256_set1_pd(Li9[i+8]);
        yi9 = sub256d(bi9,yi9);
        yi9 = fnmadd256d(Li9_8, yi8, yi9);

        __m256d bi10 = _mm256_blend_pd(bi8v, zero, 0b1011);
        __m256d Li10_8 = _mm256_set1_pd(Li10[i+8]);
        __m256d Li10_9 = _mm256_set1_pd(Li10[i+9]);
        yi10 = sub256d(bi10,yi10);
        yi10 = fnmadd256d(Li10_8, yi8, yi10);
        yi10 = fnmadd256d(Li10_9, yi9, yi10);

        __m256d bi11 = _mm256_blend_pd(bi8v, zero, 0b0111);
        __m256d Li11_8 = _mm256_set1_pd(Li11[i+8]);
        __m256d Li11_9 = _mm256_set1_pd(Li11[i+9]);
        __m256d Li11_10 = _mm256_set1_pd(Li11[i+10]);
        yi11 = sub256d(bi11,yi11);
        yi11 = fnmadd256d(Li11_8, yi8, yi11);
        yi11 = fnmadd256d(Li11_9, yi9, yi11);
        yi11 = fnmadd256d(Li11_10, yi10, yi11);

        __m256d yi8v = hsum4x4(yi8,yi9,yi10,yi11);
        store256d(&b[i+8],yi8v);
    }

    if(remain>=8)
    {
        int i=this->N-remain;
        const double* __restrict__ Li0 = &v_[(i+0)*alignN+0];
        const double* __restrict__ Li1 = &v_[(i+1)*alignN+0];
        const double* __restrict__ Li2 = &v_[(i+2)*alignN+0];
        const double* __restrict__ Li3 = &v_[(i+3)*alignN+0];
        const double* __restrict__ Li4 = &v_[(i+4)*alignN+0];
        const double* __restrict__ Li5 = &v_[(i+5)*alignN+0];
        const double* __restrict__ Li6 = &v_[(i+6)*alignN+0];
        const double* __restrict__ Li7 = &v_[(i+7)*alignN+0];
        __m256d yi0 = _mm256_setzero_pd();
        __m256d yi1 = _mm256_setzero_pd();
        __m256d yi2 = _mm256_setzero_pd();
        __m256d yi3 = _mm256_setzero_pd();
        __m256d yi4 = _mm256_setzero_pd();
        __m256d yi5 = _mm256_setzero_pd();
        __m256d yi6 = _mm256_setzero_pd();
        __m256d yi7 = _mm256_setzero_pd();
        for(int j = 0;j<i;j=j+12)
        {
            __m256d yjv = load256d(&b[j]);
            yi0 = fmadd256d(load256d(&Li0[j]),yjv,yi0);
            yi1 = fmadd256d(load256d(&Li1[j]),yjv,yi1);
            yi2 = fmadd256d(load256d(&Li2[j]),yjv,yi2);
            yi3 = fmadd256d(load256d(&Li3[j]),yjv,yi3);
            yi4 = fmadd256d(load256d(&Li4[j]),yjv,yi4);
            yi5 = fmadd256d(load256d(&Li5[j]),yjv,yi5);
            yi6 = fmadd256d(load256d(&Li6[j]),yjv,yi6);
            yi7 = fmadd256d(load256d(&Li7[j]),yjv,yi7);

            yjv = load256d(&b[j+4]);
            yi0 = fmadd256d(load256d(&Li0[j+4]),yjv,yi0);
            yi1 = fmadd256d(load256d(&Li1[j+4]),yjv,yi1);
            yi2 = fmadd256d(load256d(&Li2[j+4]),yjv,yi2);
            yi3 = fmadd256d(load256d(&Li3[j+4]),yjv,yi3);
            yi4 = fmadd256d(load256d(&Li4[j+4]),yjv,yi4);
            yi5 = fmadd256d(load256d(&Li5[j+4]),yjv,yi5);
            yi6 = fmadd256d(load256d(&Li6[j+4]),yjv,yi6);
            yi7 = fmadd256d(load256d(&Li7[j+4]),yjv,yi7);

            yjv = load256d(&b[j+8]);
            yi0 = fmadd256d(load256d(&Li0[j+8]),yjv,yi0);
            yi1 = fmadd256d(load256d(&Li1[j+8]),yjv,yi1);
            yi2 = fmadd256d(load256d(&Li2[j+8]),yjv,yi2);
            yi3 = fmadd256d(load256d(&Li3[j+8]),yjv,yi3);
            yi4 = fmadd256d(load256d(&Li4[j+8]),yjv,yi4);
            yi5 = fmadd256d(load256d(&Li5[j+8]),yjv,yi5);
            yi6 = fmadd256d(load256d(&Li6[j+8]),yjv,yi6);
            yi7 = fmadd256d(load256d(&Li7[j+8]),yjv,yi7);
        }
        __m256d zero = _mm256_setzero_pd();

        // Update 1st 4x4 diagonal block
        __m256d bi0v = load256d(&b[i+0]);

        __m256d bi0 = _mm256_blend_pd(bi0v, zero, 0b1110);
        yi0 = sub256d(bi0,yi0);

        __m256d bi1 = _mm256_blend_pd(bi0v, zero, 0b1101);
        __m256d Li10 = _mm256_set1_pd(Li1[i+0]);
        yi1 = sub256d(bi1,yi1);
        yi1 = fnmadd256d(Li10, yi0, yi1);

        __m256d bi2 = _mm256_blend_pd(bi0v, zero, 0b1011);
        __m256d Li20 = _mm256_set1_pd(Li2[i+0]);
        __m256d Li21 = _mm256_set1_pd(Li2[i+1]);
        yi2 = sub256d(bi2,yi2);
        yi2 = fnmadd256d(Li20, yi0, yi2);
        yi2 = fnmadd256d(Li21, yi1, yi2);

        __m256d bi3 = _mm256_blend_pd(bi0v, zero, 0b0111);
        __m256d Li30 = _mm256_set1_pd(Li3[i+0]);
        __m256d Li31 = _mm256_set1_pd(Li3[i+1]);
        __m256d Li32 = _mm256_set1_pd(Li3[i+2]); 
        yi3 = sub256d(bi3,yi3);
        yi3 = fnmadd256d(Li30, yi0, yi3);
        yi3 = fnmadd256d(Li31, yi1, yi3);
        yi3 = fnmadd256d(Li32, yi2, yi3);


        __m256d yi0v = hsum4x4(yi0,yi1,yi2,yi3);
        store256d(&b[i],yi0v);

        // update 2nd 4x4 block
        yi4 = fmadd256d(load256d(&Li4[i]),yi0v,yi4);
        yi5 = fmadd256d(load256d(&Li5[i]),yi0v,yi5);
        yi6 = fmadd256d(load256d(&Li6[i]),yi0v,yi6);
        yi7 = fmadd256d(load256d(&Li7[i]),yi0v,yi7);

        // update 2nd 4x4 diagonal block
        __m256d bi4v = load256d(&b[i+4]);

        __m256d bi4 = _mm256_blend_pd(bi4v, zero, 0b1110);
        yi4 = sub256d(bi4,yi4);

        __m256d bi5 = _mm256_blend_pd(bi4v, zero, 0b1101);
        __m256d Li54 = _mm256_set1_pd(Li5[i+4]);
        yi5 = sub256d(bi5,yi5);
        yi5 = fnmadd256d(Li54, yi4, yi5); 

        __m256d bi6 = _mm256_blend_pd(bi4v, zero, 0b1011);
        __m256d Li64 = _mm256_set1_pd(Li6[i+4]); 
        __m256d Li65 = _mm256_set1_pd(Li6[i+5]);
        yi6 = sub256d(bi6,yi6);
        yi6 = fnmadd256d(Li64, yi4, yi6); 
        yi6 = fnmadd256d(Li65, yi5, yi6);

        __m256d bi7 = _mm256_blend_pd(bi4v, zero, 0b0111);
        __m256d Li74 = _mm256_set1_pd(Li7[i+4]);
        __m256d Li75 = _mm256_set1_pd(Li7[i+5]);
        __m256d Li76 = _mm256_set1_pd(Li7[i+6]);
        yi7 = sub256d(bi7,yi7);
        yi7 = fnmadd256d(Li74, yi4, yi7);
        yi7 = fnmadd256d(Li75, yi5, yi7); 
        yi7 = fnmadd256d(Li76, yi6, yi7);


        __m256d yi4v = hsum4x4(yi4,yi5,yi6,yi7); 
        store256d(&b[i+4],yi4v);
    }

    if(remain>=4&&remain<8)
    {
        int i=this->N-remain;
        const double* __restrict__ Li0 = &v_[(i+0)*alignN+0];
        const double* __restrict__ Li1 = &v_[(i+1)*alignN+0];
        const double* __restrict__ Li2 = &v_[(i+2)*alignN+0];
        const double* __restrict__ Li3 = &v_[(i+3)*alignN+0];
        __m256d yi0 = _mm256_setzero_pd();
        __m256d yi1 = _mm256_setzero_pd();
        __m256d yi2 = _mm256_setzero_pd();
        __m256d yi3 = _mm256_setzero_pd();
        for(int j = 0;j<i;j=j+12)
        {
            __m256d yjv = load256d(&b[j]);
            yi0 = fmadd256d(load256d(&Li0[j]),yjv,yi0);
            yi1 = fmadd256d(load256d(&Li1[j]),yjv,yi1);
            yi2 = fmadd256d(load256d(&Li2[j]),yjv,yi2);
            yi3 = fmadd256d(load256d(&Li3[j]),yjv,yi3);

            yjv = load256d(&b[j+4]);
            yi0 = fmadd256d(load256d(&Li0[j+4]),yjv,yi0);
            yi1 = fmadd256d(load256d(&Li1[j+4]),yjv,yi1);
            yi2 = fmadd256d(load256d(&Li2[j+4]),yjv,yi2);
            yi3 = fmadd256d(load256d(&Li3[j+4]),yjv,yi3);

            yjv = load256d(&b[j+8]);
            yi0 = fmadd256d(load256d(&Li0[j+8]),yjv,yi0);
            yi1 = fmadd256d(load256d(&Li1[j+8]),yjv,yi1);
            yi2 = fmadd256d(load256d(&Li2[j+8]),yjv,yi2);
            yi3 = fmadd256d(load256d(&Li3[j+8]),yjv,yi3);
        }
        __m256d zero = _mm256_setzero_pd();

        // Update the 4x4 diagonal block
        __m256d bi0v = load256d(&b[i+0]);

        __m256d bi0 = _mm256_blend_pd(bi0v, zero, 0b1110);
        yi0 = sub256d(bi0,yi0);

        __m256d bi1 = _mm256_blend_pd(bi0v, zero, 0b1101);
        __m256d Li10 = _mm256_set1_pd(Li1[i+0]);
        yi1 = sub256d(bi1,yi1);
        yi1 = fnmadd256d(Li10, yi0, yi1); 

        __m256d bi2 = _mm256_blend_pd(bi0v, zero, 0b1011); 
        __m256d Li20 = _mm256_set1_pd(Li2[i+0]);
        __m256d Li21 = _mm256_set1_pd(Li2[i+1]);
        yi2 = sub256d(bi2,yi2);
        yi2 = fnmadd256d(Li20, yi0, yi2);
        yi2 = fnmadd256d(Li21, yi1, yi2);

        __m256d bi3 = _mm256_blend_pd(bi0v, zero, 0b0111);
        __m256d Li30 = _mm256_set1_pd(Li3[i+0]);
        __m256d Li31 = _mm256_set1_pd(Li3[i+1]);
        __m256d Li32 = _mm256_set1_pd(Li3[i+2]);
        yi3 = sub256d(bi3,yi3);
        yi3 = fnmadd256d(Li30, yi0, yi3);
        yi3 = fnmadd256d(Li31, yi1, yi3);
        yi3 = fnmadd256d(Li32, yi2, yi3);


        __m256d yi0v = hsum4x4(yi0,yi1,yi2,yi3);
        store256d(&b[i],yi0v);
    }

    {


        int nRows = (remain)%4;

        for(int i=N-nRows;i<N;i++)
        {
            const double* __restrict__ L = &v_[(i+0)*alignN+0];

            __m256d yv = _mm256_setzero_pd();
            double ys = 0;

            //nRowsx12 block from main loop
            for(int j = 0;j<N-remain;j=j+12)
            {
                __m256d
                yjv = load256d(&b[j+0]);
                yv = fmadd256d(load256d(&L[j]),yjv,yv);

                yjv = load256d(&b[j+4]);
                yv = fmadd256d(load256d(&L[j+4]),yjv,yv);

                yjv = load256d(&b[j+8]);
                yv = fmadd256d(load256d(&L[j+8]),yjv,yv);
            }
            for (int j = N - remain; j < N - nRows; j += 4)
            {
                __m256d yjv = load256d(&b[j]);
                yv = fmadd256d(load256d(&L[j]), yjv, yv);
            }
            for(int j = N-nRows;j<i;j++)
            {
                ys = ys + L[j]*b[j];
            }
            b[i] = b[i] - hsum4(yv) - ys;
        }
    }

    int tail = N%12;

    for(int i = N-12; i >= 0; i=i-12)
    {
        const double* __restrict__ U11 = &this->v_[(i+11)*alignN];
        const double* __restrict__ U10 = &this->v_[(i+10)*alignN];
        const double* __restrict__ U9 = &this->v_[(i+9)*alignN];
        const double* __restrict__ U8 = &this->v_[(i+8)*alignN];
        const double* __restrict__ U7 = &this->v_[(i+7)*alignN];
        const double* __restrict__ U6 = &this->v_[(i+6)*alignN];
        const double* __restrict__ U5 = &this->v_[(i+5)*alignN];
        const double* __restrict__ U4 = &this->v_[(i+4)*alignN];
        const double* __restrict__ U3 = &this->v_[(i+3)*alignN];
        const double* __restrict__ U2 = &this->v_[(i+2)*alignN];
        const double* __restrict__ U1 = &this->v_[(i+1)*alignN];
        const double* __restrict__ U0 = &this->v_[(i+0)*alignN];

        __m256d x11 = _mm256_setzero_pd();
        __m256d x10 = _mm256_setzero_pd();
        __m256d x9 = _mm256_setzero_pd();
        __m256d x8 = _mm256_setzero_pd();
        __m256d x7 = _mm256_setzero_pd();
        __m256d x6 = _mm256_setzero_pd();
        __m256d x5 = _mm256_setzero_pd();
        __m256d x4 = _mm256_setzero_pd();
        __m256d x3 = _mm256_setzero_pd();
        __m256d x2 = _mm256_setzero_pd();
        __m256d x1 = _mm256_setzero_pd();
        __m256d x0 = _mm256_setzero_pd();

        for(int j = i+12; j <= N-12; j=j+12)
        {
            __m256d x12v = loadu256d(&b[j+0]);
            x11 = fmadd256d(loadu256d(&U11[j+0]),x12v,x11);
            x10 = fmadd256d(loadu256d(&U10[j+0]),x12v,x10);
            x9 = fmadd256d(loadu256d(&U9[j+0]),x12v,x9);
            x8 = fmadd256d(loadu256d(&U8[j+0]),x12v,x8);
            x7 = fmadd256d(loadu256d(&U7[j+0]),x12v,x7);
            x6 = fmadd256d(loadu256d(&U6[j+0]),x12v,x6);
            x5 = fmadd256d(loadu256d(&U5[j+0]),x12v,x5);
            x4 = fmadd256d(loadu256d(&U4[j+0]),x12v,x4);
            x3 = fmadd256d(loadu256d(&U3[j+0]),x12v,x3);
            x2 = fmadd256d(loadu256d(&U2[j+0]),x12v,x2);
            x1 = fmadd256d(loadu256d(&U1[j+0]),x12v,x1);
            x0 = fmadd256d(loadu256d(&U0[j+0]),x12v,x0);

            __m256d x16v = loadu256d(&b[j+4]);
            x11 = fmadd256d(loadu256d(&U11[j+4]),x16v,x11);
            x10 = fmadd256d(loadu256d(&U10[j+4]),x16v,x10);
            x9 = fmadd256d(loadu256d(&U9[j+4]),x16v,x9);
            x8 = fmadd256d(loadu256d(&U8[j+4]),x16v,x8);
            x7 = fmadd256d(loadu256d(&U7[j+4]),x16v,x7);
            x6 = fmadd256d(loadu256d(&U6[j+4]),x16v,x6);
            x5 = fmadd256d(loadu256d(&U5[j+4]),x16v,x5);
            x4 = fmadd256d(loadu256d(&U4[j+4]),x16v,x4);
            x3 = fmadd256d(loadu256d(&U3[j+4]),x16v,x3);
            x2 = fmadd256d(loadu256d(&U2[j+4]),x16v,x2);
            x1 = fmadd256d(loadu256d(&U1[j+4]),x16v,x1);
            x0 = fmadd256d(loadu256d(&U0[j+4]),x16v,x0);

            __m256d x20v = loadu256d(&b[j+8]);
            x11 = fmadd256d(loadu256d(&U11[j+8]),x20v,x11);
            x10 = fmadd256d(loadu256d(&U10[j+8]),x20v,x10);
            x9 = fmadd256d(loadu256d(&U9[j+8]),x20v,x9);
            x8 = fmadd256d(loadu256d(&U8[j+8]),x20v,x8);
            x7 = fmadd256d(loadu256d(&U7[j+8]),x20v,x7);
            x6 = fmadd256d(loadu256d(&U6[j+8]),x20v,x6);
            x5 = fmadd256d(loadu256d(&U5[j+8]),x20v,x5);
            x4 = fmadd256d(loadu256d(&U4[j+8]),x20v,x4);
            x3 = fmadd256d(loadu256d(&U3[j+8]),x20v,x3);
            x2 = fmadd256d(loadu256d(&U2[j+8]),x20v,x2);
            x1 = fmadd256d(loadu256d(&U1[j+8]),x20v,x1);
            x0 = fmadd256d(loadu256d(&U0[j+8]),x20v,x0);
        }

        __m256d zero = _mm256_setzero_pd();
        __m256d b8v = loadu256d(&b[i+8]);

        __m256d b11 = _mm256_blend_pd(b8v, zero, 0b0111);
        x11 = sub256d(b11,x11);
        __m256d ru8v = loadu256d(&this->invD[i+8]);
        __m256d ru11 = b11;
        ru11 = _mm256_permute4x64_pd(ru8v,0b11111111);
        x11 = mul256d(x11,ru11);

        __m256d U10_11 = _mm256_set1_pd(U10[i+11]);
        x10 = fmadd256d(U10_11,x11,x10);
        __m256d b10 = _mm256_blend_pd(b8v, zero, 0b1011);
        x10 = sub256d(b10,x10);
        __m256d ru10 = _mm256_permute4x64_pd(ru8v,0b10101010);
        x10 = mul256d(x10,ru10);

        __m256d U9_10 = _mm256_set1_pd(U9[i+10]);
        x9 = fmadd256d(U9_10,x10,x9);
        __m256d U9_11 = _mm256_set1_pd(U9[i+11]);
        x9 = fmadd256d(U9_11,x11,x9);
        __m256d b9 = _mm256_blend_pd(b8v, zero, 0b1101);
        x9 = sub256d(b9,x9);
        __m256d ru9 = _mm256_permute4x64_pd(ru8v,0b01010101);
        x9 = mul256d(x9,ru9);

        __m256d U89 = _mm256_set1_pd(U8[i+9]);
        x8 = fmadd256d(U89,x9,x8);
        __m256d U810 = _mm256_set1_pd(U8[i+10]);
        x8 = fmadd256d(U810,x10,x8);
        __m256d U811 = _mm256_set1_pd(U8[i+11]);
        x8 = fmadd256d(U811,x11,x8);
        __m256d b8 = _mm256_blend_pd(b8v, zero, 0b1110);
        x8 = sub256d(b8,x8);
        __m256d ru8 = _mm256_permute4x64_pd(ru8v,0b00000000);
        x8 = mul256d(x8,ru8);

        x8 = _mm256_hadd_pd(x8, x9);
        x9 = _mm256_hadd_pd(x10, x11);
        x8 = _mm256_permute4x64_pd (x8, _MM_SHUFFLE(3,1,2,0));
        x9 = _mm256_permute4x64_pd(x9, _MM_SHUFFLE(3,1,2,0));
        x8 = _mm256_hadd_pd(x8, x9);
        __m256d x8v = _mm256_permute4x64_pd(x8, _MM_SHUFFLE(3,1,2,0));
        storeu256d(&b[i+8],x8v);

        x7 = fmadd256d(loadu256d(&U7[i+8]),x8v,x7);
        x6 = fmadd256d(loadu256d(&U6[i+8]),x8v,x6);
        x5 = fmadd256d(loadu256d(&U5[i+8]),x8v,x5);
        x4 = fmadd256d(loadu256d(&U4[i+8]),x8v,x4);
        x3 = fmadd256d(loadu256d(&U3[i+8]),x8v,x3);
        x2 = fmadd256d(loadu256d(&U2[i+8]),x8v,x2);
        x1 = fmadd256d(loadu256d(&U1[i+8]),x8v,x1);
        x0 = fmadd256d(loadu256d(&U0[i+8]),x8v,x0);

        __m256d b4v = loadu256d(&b[i+4]);
        __m256d ru4v = loadu256d(&this->invD[i+4]);

        __m256d b7 = _mm256_blend_pd(b4v, zero, 0b0111);
        __m256d ru7 = _mm256_permute4x64_pd(ru4v,0b11111111);
        x7 = sub256d(b7,x7);
        x7 = mul256d(x7,ru7);

        __m256d U67 = _mm256_set1_pd(U6[i+7]);
        x6 = fmadd256d(U67,x7,x6);
        __m256d b6 = _mm256_blend_pd(b4v, zero, 0b1011);
        __m256d ru6 = _mm256_permute4x64_pd(ru4v,0b10101010);
        x6 = sub256d(b6,x6);
        x6 = mul256d(x6,ru6);

        __m256d U56 = _mm256_set1_pd(U5[i+6]);
        __m256d U57 = _mm256_set1_pd(U5[i+7]);
        x5 = fmadd256d(U56,x6,x5);
        x5 = fmadd256d(U57,x7,x5);
        __m256d b5 = _mm256_blend_pd(b4v, zero, 0b1101);
        __m256d ru5 = _mm256_permute4x64_pd(ru4v,0b01010101);
        x5 = sub256d(b5,x5);
        x5 = mul256d(x5,ru5);

        __m256d U45 = _mm256_set1_pd(U4[i+5]);
        __m256d U46 = _mm256_set1_pd(U4[i+6]);
        __m256d U47 = _mm256_set1_pd(U4[i+7]);
        x4 = fmadd256d(U45,x5,x4);
        x4 = fmadd256d(U46,x6,x4);
        x4 = fmadd256d(U47,x7,x4);
        __m256d b4 = _mm256_blend_pd(b4v, zero, 0b1110);
        __m256d ru4 = _mm256_permute4x64_pd(ru4v,0b00000000);
        x4 = sub256d(b4,x4);
        x4 = mul256d(x4,ru4);

        __m256d x4v = hsum4x4(x4,x5,x6,x7);
        storeu256d(&b[i+4],x4v);

        // update 3rd 4x4 block
        x3 = fmadd256d(loadu256d(&U3[i+4]),x4v,x3);
        x2 = fmadd256d(loadu256d(&U2[i+4]),x4v,x2);
        x1 = fmadd256d(loadu256d(&U1[i+4]),x4v,x1);
        x0 = fmadd256d(loadu256d(&U0[i+4]),x4v,x0);

        __m256d b0v = loadu256d(&b[i+0]);
        __m256d ru0v = loadu256d(&this->invD[i+0]);

        __m256d b3 = _mm256_blend_pd(b0v, zero, 0b0111);
        __m256d ru3 = _mm256_permute4x64_pd(ru0v,0b11111111);
        x3 = sub256d(b3,x3);
        x3 = mul256d(x3,ru3);

        __m256d U23 = _mm256_set1_pd(U2[i+3]);
        x2 = fmadd256d(U23,x3,x2);
        __m256d b2 = _mm256_blend_pd(b0v, zero, 0b1011);
        __m256d ru2 = _mm256_permute4x64_pd(ru0v,0b10101010);
        x2 = sub256d(b2,x2);
        x2 = mul256d(x2,ru2);

        __m256d U12 = _mm256_set1_pd(U1[i+2]);
        __m256d U13 = _mm256_set1_pd(U1[i+3]);
        x1 = fmadd256d(U12,x2,x1);
        x1 = fmadd256d(U13,x3,x1);
        __m256d b1 = _mm256_blend_pd(b0v, zero, 0b1101);
        __m256d ru1 = _mm256_permute4x64_pd(ru0v,0b01010101);
        x1 = sub256d(b1,x1);
        x1 = mul256d(x1,ru1);

        __m256d U01 = _mm256_set1_pd(U0[i+1]);
        __m256d U02 = _mm256_set1_pd(U0[i+2]);
        __m256d U03 = _mm256_set1_pd(U0[i+3]);
        x0 = fmadd256d(U01,x1,x0);
        x0 = fmadd256d(U02,x2,x0);
        x0 = fmadd256d(U03,x3,x0);
        __m256d b0 = _mm256_blend_pd(b0v, zero, 0b1110);
        __m256d ru0 = _mm256_permute4x64_pd(ru0v,0b00000000);
        x0 = sub256d(b0,x0);
        x0 = mul256d(x0,ru0);

        __m256d x0v = hsum4x4(x0,x1,x2,x3);
        storeu256d(&b[i+0],x0v);
    }

    // tail=8,9,10,11. Use 8x8 block
    if(tail>=8)
    {
        int i = tail - 8;

        const double* __restrict__ U7 = &this->v_[(i+7)*alignN];
        const double* __restrict__ U6 = &this->v_[(i+6)*alignN];
        const double* __restrict__ U5 = &this->v_[(i+5)*alignN];
        const double* __restrict__ U4 = &this->v_[(i+4)*alignN];
        const double* __restrict__ U3 = &this->v_[(i+3)*alignN];
        const double* __restrict__ U2 = &this->v_[(i+2)*alignN];
        const double* __restrict__ U1 = &this->v_[(i+1)*alignN];
        const double* __restrict__ U0 = &this->v_[(i+0)*alignN];

        __m256d x7 = _mm256_setzero_pd();
        __m256d x6 = _mm256_setzero_pd();
        __m256d x5 = _mm256_setzero_pd();
        __m256d x4 = _mm256_setzero_pd();
        __m256d x3 = _mm256_setzero_pd();
        __m256d x2 = _mm256_setzero_pd();
        __m256d x1 = _mm256_setzero_pd();
        __m256d x0 = _mm256_setzero_pd();

        int NN = N;
        for(int j = i+8; j <= NN-12; j=j+12)
        {

            __m256d x12v = loadu256d(&b[j+0]);
            x7 = fmadd256d(loadu256d(&U7[j+0]),x12v,x7);
            x6 = fmadd256d(loadu256d(&U6[j+0]),x12v,x6);
            x5 = fmadd256d(loadu256d(&U5[j+0]),x12v,x5);
            x4 = fmadd256d(loadu256d(&U4[j+0]),x12v,x4);
            x3 = fmadd256d(loadu256d(&U3[j+0]),x12v,x3);
            x2 = fmadd256d(loadu256d(&U2[j+0]),x12v,x2);
            x1 = fmadd256d(loadu256d(&U1[j+0]),x12v,x1);
            x0 = fmadd256d(loadu256d(&U0[j+0]),x12v,x0);

            __m256d x16v = loadu256d(&b[j+4]);
            x7 = fmadd256d(loadu256d(&U7[j+4]),x16v,x7);
            x6 = fmadd256d(loadu256d(&U6[j+4]),x16v,x6);
            x5 = fmadd256d(loadu256d(&U5[j+4]),x16v,x5);
            x4 = fmadd256d(loadu256d(&U4[j+4]),x16v,x4);
            x3 = fmadd256d(loadu256d(&U3[j+4]),x16v,x3);
            x2 = fmadd256d(loadu256d(&U2[j+4]),x16v,x2);
            x1 = fmadd256d(loadu256d(&U1[j+4]),x16v,x1);
            x0 = fmadd256d(loadu256d(&U0[j+4]),x16v,x0);

            __m256d x20v = loadu256d(&b[j+8]);
            x7 = fmadd256d(loadu256d(&U7[j+8]),x20v,x7);
            x6 = fmadd256d(loadu256d(&U6[j+8]),x20v,x6);
            x5 = fmadd256d(loadu256d(&U5[j+8]),x20v,x5);
            x4 = fmadd256d(loadu256d(&U4[j+8]),x20v,x4);
            x3 = fmadd256d(loadu256d(&U3[j+8]),x20v,x3);
            x2 = fmadd256d(loadu256d(&U2[j+8]),x20v,x2);
            x1 = fmadd256d(loadu256d(&U1[j+8]),x20v,x1);
            x0 = fmadd256d(loadu256d(&U0[j+8]),x20v,x0);
        }
        __m256d zero = _mm256_setzero_pd();
        // update 1st 4x4 block
        __m256d b4v = loadu256d(&b[i+4]);
        __m256d ru4v = loadu256d(&this->invD[i+4]);

        __m256d b7 = _mm256_blend_pd(b4v, zero, 0b0111);
        __m256d ru7 = _mm256_permute4x64_pd(ru4v,0b11111111);
        x7 = sub256d(b7,x7);
        x7 = mul256d(x7,ru7);

        __m256d U67 = _mm256_set1_pd(U6[i+7]);
        x6 = fmadd256d(U67,x7,x6);
        __m256d b6 = _mm256_blend_pd(b4v, zero, 0b1011);
        __m256d ru6 = _mm256_permute4x64_pd(ru4v,0b10101010);
        x6 = sub256d(b6,x6);
        x6 = mul256d(x6,ru6);

        __m256d U56 = _mm256_set1_pd(U5[i+6]);
        __m256d U57 = _mm256_set1_pd(U5[i+7]);
        x5 = fmadd256d(U56,x6,x5);
        x5 = fmadd256d(U57,x7,x5);
        __m256d b5 = _mm256_blend_pd(b4v, zero, 0b1101);
        __m256d ru5 = _mm256_permute4x64_pd(ru4v,0b01010101);
        x5 = sub256d(b5,x5);
        x5 = mul256d(x5,ru5);

        __m256d U45 = _mm256_set1_pd(U4[i+5]);
        __m256d U46 = _mm256_set1_pd(U4[i+6]);
        __m256d U47 = _mm256_set1_pd(U4[i+7]);
        x4 = fmadd256d(U45,x5,x4);
        x4 = fmadd256d(U46,x6,x4);
        x4 = fmadd256d(U47,x7,x4);
        __m256d b4 = _mm256_blend_pd(b4v, zero, 0b1110);
        __m256d ru4 = _mm256_permute4x64_pd(ru4v,0b00000000);
        x4 = sub256d(b4,x4);
        x4 = mul256d(x4,ru4);

        __m256d x4v = hsum4x4(x4,x5,x6,x7);
        storeu256d(&b[i+4],x4v);

        // update 3rd 4x4 block
        x3 = fmadd256d(loadu256d(&U3[i+4]),x4v,x3);
        x2 = fmadd256d(loadu256d(&U2[i+4]),x4v,x2);
        x1 = fmadd256d(loadu256d(&U1[i+4]),x4v,x1);
        x0 = fmadd256d(loadu256d(&U0[i+4]),x4v,x0);

        __m256d b0v = loadu256d(&b[i+0]);
        __m256d ru0v = loadu256d(&this->invD[i+0]);

        __m256d b3 = _mm256_blend_pd(b0v, zero, 0b0111);  //[0,0,0,b3]
        __m256d ru3 = _mm256_permute4x64_pd(ru0v,0b11111111);
        x3 = sub256d(b3,x3);
        x3 = mul256d(x3,ru3);

        __m256d U23 = _mm256_set1_pd(U2[i+3]);
        x2 = fmadd256d(U23,x3,x2);
        __m256d b2 = _mm256_blend_pd(b0v, zero, 0b1011);  //[0,b2,0,0]
        __m256d ru2 = _mm256_permute4x64_pd(ru0v,0b10101010);
        x2 = sub256d(b2,x2);
        x2 = mul256d(x2,ru2);

        __m256d U12 = _mm256_set1_pd(U1[i+2]);
        __m256d U13 = _mm256_set1_pd(U1[i+3]);
        x1 = fmadd256d(U12,x2,x1);
        x1 = fmadd256d(U13,x3,x1);
        __m256d b1 = _mm256_blend_pd(b0v, zero, 0b1101);  //[0,b1,0,0]
        __m256d ru1 = _mm256_permute4x64_pd(ru0v,0b01010101);
        x1 = sub256d(b1,x1);
        x1 = mul256d(x1,ru1);

        __m256d U01 = _mm256_set1_pd(U0[i+1]);
        __m256d U02 = _mm256_set1_pd(U0[i+2]);
        __m256d U03 = _mm256_set1_pd(U0[i+3]);
        x0 = fmadd256d(U01,x1,x0);
        x0 = fmadd256d(U02,x2,x0);
        x0 = fmadd256d(U03,x3,x0);
        __m256d b0 = _mm256_blend_pd(b0v, zero, 0b1110);  //[b0,0,0,0]
        __m256d ru0 = _mm256_permute4x64_pd(ru0v,0b00000000);
        x0 = sub256d(b0,x0);
        x0 = mul256d(x0,ru0);

        __m256d x0v = hsum4x4(x0,x1,x2,x3);
        storeu256d(&b[i+0],x0v);
    }

    else if(tail>=4 && tail<8)
    {
        int i = tail - 4;
        const double* __restrict__ U3 = &this->v_[(i+3)*alignN];
        const double* __restrict__ U2 = &this->v_[(i+2)*alignN];
        const double* __restrict__ U1 = &this->v_[(i+1)*alignN];
        const double* __restrict__ U0 = &this->v_[(i+0)*alignN];

        __m256d x3 = _mm256_setzero_pd();
        __m256d x2 = _mm256_setzero_pd();
        __m256d x1 = _mm256_setzero_pd();
        __m256d x0 = _mm256_setzero_pd();

        for(int j = i+4; j <= N-12; j=j+12)
        {
            __m256d x12v = loadu256d(&b[j+0]);
            x3 = fmadd256d(loadu256d(&U3[j+0]),x12v,x3);
            x2 = fmadd256d(loadu256d(&U2[j+0]),x12v,x2);
            x1 = fmadd256d(loadu256d(&U1[j+0]),x12v,x1);
            x0 = fmadd256d(loadu256d(&U0[j+0]),x12v,x0);

            __m256d x16v = loadu256d(&b[j+4]);
            x3 = fmadd256d(loadu256d(&U3[j+4]),x16v,x3);
            x2 = fmadd256d(loadu256d(&U2[j+4]),x16v,x2);
            x1 = fmadd256d(loadu256d(&U1[j+4]),x16v,x1);
            x0 = fmadd256d(loadu256d(&U0[j+4]),x16v,x0);

            __m256d x20v = loadu256d(&b[j+8]);
            x3 = fmadd256d(loadu256d(&U3[j+8]),x20v,x3);
            x2 = fmadd256d(loadu256d(&U2[j+8]),x20v,x2);
            x1 = fmadd256d(loadu256d(&U1[j+8]),x20v,x1);
            x0 = fmadd256d(loadu256d(&U0[j+8]),x20v,x0);
        }
        __m256d zero = _mm256_setzero_pd();
        __m256d b0v = loadu256d(&b[i+0]);
        __m256d ru0v = loadu256d(&this->invD[i+0]);

        __m256d b3 = _mm256_blend_pd(b0v, zero, 0b0111);  //[0,0,0,b3]
        __m256d ru3 = _mm256_permute4x64_pd(ru0v,0b11111111);
        x3 = sub256d(b3,x3);
        x3 = mul256d(x3,ru3);

        __m256d U23 = _mm256_set1_pd(U2[i+3]);
        x2 = fmadd256d(U23,x3,x2);
        __m256d b2 = _mm256_blend_pd(b0v, zero, 0b1011);  //[0,b2,0,0]
        __m256d ru2 = _mm256_permute4x64_pd(ru0v,0b10101010);
        x2 = sub256d(b2,x2);
        x2 = mul256d(x2,ru2);

        __m256d U12 = _mm256_set1_pd(U1[i+2]);
        __m256d U13 = _mm256_set1_pd(U1[i+3]);
        x1 = fmadd256d(U12,x2,x1);
        x1 = fmadd256d(U13,x3,x1);
        __m256d b1 = _mm256_blend_pd(b0v, zero, 0b1101);  //[0,b1,0,0]
        __m256d ru1 = _mm256_permute4x64_pd(ru0v,0b01010101);
        x1 = sub256d(b1,x1);
        x1 = mul256d(x1,ru1);

        __m256d U01 = _mm256_set1_pd(U0[i+1]);
        __m256d U02 = _mm256_set1_pd(U0[i+2]);
        __m256d U03 = _mm256_set1_pd(U0[i+3]);
        x0 = fmadd256d(U01,x1,x0);
        x0 = fmadd256d(U02,x2,x0);
        x0 = fmadd256d(U03,x3,x0);
        __m256d b0 = _mm256_blend_pd(b0v, zero, 0b1110);  //[b0,0,0,0]
        __m256d ru0 = _mm256_permute4x64_pd(ru0v,0b00000000);
        x0 = sub256d(b0,x0);
        x0 = mul256d(x0,ru0);

        __m256d x0v = hsum4x4(x0,x1,x2,x3);
        storeu256d(&b[i+0],x0v);
    }

    // row by row if tail has 1 or 2 or 3 rows left after 8x8 and 4x4 block
    {
        const int k =((tail >> 2) << 2) + 1;
        const int istart = tail - k;
        const int jstart = tail;
        for(int i=istart;i>=0;i--)
        {

            const double* __restrict__ U = &this->v_[(i)*alignN];
            __m256d x = _mm256_setzero_pd();
            double xs = 0;
            int NN = N;
            for(int j = jstart; j <= NN-12; j=j+12)
            {
                __m256d x12v = loadu256d(&b[j+0]);
                x = fmadd256d(loadu256d(&U[j+0]),x12v,x);
                __m256d x16v = loadu256d(&b[j+4]);
                x = fmadd256d(loadu256d(&U[j+4]),x16v,x);
                __m256d x20v = loadu256d(&b[j+8]);
                x = fmadd256d(loadu256d(&U[j+8]),x20v,x);
            }
            if(tail>=8)
            {
                int j = tail - 8;
                __m256d x12v = loadu256d(&b[j+0]);
                x = fmadd256d(loadu256d(&U[j+0]),x12v,x);
                __m256d x16v = loadu256d(&b[j+4]);
                x = fmadd256d(loadu256d(&U[j+4]),x16v,x);
            }
            else if(tail >=4 && tail < 8)
            {
                int j = tail - 4;
                __m256d x12v = loadu256d(&b[j+0]);
                x = fmadd256d(loadu256d(&U[j+0]),x12v,x);
            }
            for(int j =istart;j>i;j--)
            {
                xs = xs + U[j]*b[j];
            }

            double result = (b[i] - xs - hsum4(x))*this->invD[i];
            b[i] = result;
        }
    }
}

// Solve LUx = b (serial forward/backward substitution).
void FastChemistry::LUsolver::xSolve_Serial
(
    double* __restrict__ b
)
{
    for(int j = 0; j<this->N;j++)
    {
        if(this->pivotIndice_[j]==j)
        {
            continue;
        }
        int jTarget = this->pivotIndice_[j];
        double temp = b[jTarget];
        b[jTarget] = b[j];
        b[j] = temp;
    }

    for(int i = 0; i < this->N; i++)
    {
        double sum = b[i];
        for(int j = 0; j < i; j++)
        {
            sum = sum - v_[i*alignN+j]*b[j];
        }
        b[i] = sum;
    }

    for(int i = N-1; i >= 0; i--)
    {
        double sum = b[i];
        for(int j = i+1; j<this->N; j++)
        {
            const double xj = b[j];
            sum = sum - v_[i*alignN+j]*xj;
        }
        b[i] = sum/v_[i*alignN+i];
    }

}

// Row equilibration: divide each row by its maximum absolute value and store
// the reciprocal factors in rowScales.  The four Remain variants only differ
// in how the last 1..3 rows are handled.
void FastChemistry::LUsolver::scaleA3
(
)
{
    int remain8 = N%8;
    for(int i = 0; i < N-3; i=i+4)
    {
        double* ptr0 = &this->v_[(i+0)*this->alignN];
        double* ptr1 = &this->v_[(i+1)*this->alignN];
        double* ptr2 = &this->v_[(i+2)*this->alignN];
        double* ptr3 = &this->v_[(i+3)*this->alignN];

        __m256d maxVal0 = _mm256_setzero_pd();
        __m256d maxVal1 = _mm256_setzero_pd();
        __m256d maxVal2 = _mm256_setzero_pd();
        __m256d maxVal3 = _mm256_setzero_pd();
        const __m256d absmask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
        for(int j = 0; j < this->N-remain8; j=j+8)
        {

            __m256d absVal00 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            __m256d absVal04 = _mm256_and_pd(load256d(&ptr0[j+4]), absmask);
            maxVal0 = _mm256_max_pd(absVal00,maxVal0);
            maxVal0 = _mm256_max_pd(absVal04,maxVal0);
            __m256d absVal10 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            __m256d absVal14 = _mm256_and_pd(load256d(&ptr1[j+4]), absmask);
            maxVal1 = _mm256_max_pd(absVal10,maxVal1);
            maxVal1 = _mm256_max_pd(absVal14,maxVal1);
            __m256d absVal20 = _mm256_and_pd(load256d(&ptr2[j+0]), absmask);
            __m256d absVal24 = _mm256_and_pd(load256d(&ptr2[j+4]), absmask);
            maxVal2 = _mm256_max_pd(absVal20,maxVal2);
            maxVal2 = _mm256_max_pd(absVal24,maxVal2);
            __m256d absVal30 = _mm256_and_pd(load256d(&ptr3[j+0]), absmask);
            __m256d absVal34 = _mm256_and_pd(load256d(&ptr3[j+4]), absmask);
            maxVal3 = _mm256_max_pd(absVal30,maxVal3);
            maxVal3 = _mm256_max_pd(absVal34,maxVal3);
        }
        if(remain8==7)
        {
            int j = this->N-7;
            __m256d absVal0 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
            __m256d absVal2 = _mm256_and_pd(load256d(&ptr2[j+0]), absmask);
            maxVal2 = _mm256_max_pd(absVal2,maxVal2);
            __m256d absVal3 = _mm256_and_pd(load256d(&ptr3[j+0]), absmask);
            maxVal3 = _mm256_max_pd(absVal3,maxVal3);
        }
        {
            int j = this->N-3;
            __m256d absVal0 = _mm256_and_pd(load256d012(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d012(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
            __m256d absVal2 = _mm256_and_pd(load256d012(&ptr2[j+0]), absmask);
            maxVal2 = _mm256_max_pd(absVal2,maxVal2);
            __m256d absVal3 = _mm256_and_pd(load256d012(&ptr3[j+0]), absmask);
            maxVal3 = _mm256_max_pd(absVal3,maxVal3);
        }

        maxVal0 = broadcastMax(maxVal0);
        maxVal1 = broadcastMax(maxVal1);
        maxVal2 = broadcastMax(maxVal2);
        maxVal3 = broadcastMax(maxVal3);

        __m256d temp1 = _mm256_unpacklo_pd(maxVal0, maxVal1);
        __m256d temp2 = _mm256_unpacklo_pd(maxVal2, maxVal3);
        __m256d result = _mm256_permute2f128_pd(temp1, temp2, 0x20);
        __m256d ones = _mm256_set1_pd(1.0);
        result =  _mm256_div_pd(ones, result);
        storeu256d(&this->rowScales[i+0],result);
        __m256d recipMaxVal0 = _mm256_permute4x64_pd(result, 0b00000000);
        __m256d recipMaxVal1 = _mm256_permute4x64_pd(result, 0b01010101);
        __m256d recipMaxVal2 = _mm256_permute4x64_pd(result, 0b10101010);
        __m256d recipMaxVal3 = _mm256_permute4x64_pd(result, 0b11111111);

        for(int j = 0; j < this->N-remain8; j=j+8)
        {
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d(&ptr3[j+0]),recipMaxVal3);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
            store256d(&ptr2[j+0],r2);
            store256d(&ptr3[j+0],r3);
            __m256d r4 = mul256d(load256d(&ptr0[j+4]),recipMaxVal0);
            __m256d r5 = mul256d(load256d(&ptr1[j+4]),recipMaxVal1);
            __m256d r6 = mul256d(load256d(&ptr2[j+4]),recipMaxVal2);
            __m256d r7 = mul256d(load256d(&ptr3[j+4]),recipMaxVal3);
            store256d(&ptr0[j+4],r4);
            store256d(&ptr1[j+4],r5);
            store256d(&ptr2[j+4],r6);
            store256d(&ptr3[j+4],r7);
        }
        if(remain8==7)
        {
            int j = this->N-7;
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d(&ptr3[j+0]),recipMaxVal3);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
            store256d(&ptr2[j+0],r2);
            store256d(&ptr3[j+0],r3);
        }
        {
            int j = this->N-3;
            __m256d r0 = mul256d(load256d012(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d012(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d012(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d012(&ptr3[j+0]),recipMaxVal3);
            store256d012(&ptr0[j+0],r0);
            store256d012(&ptr1[j+0],r1);
            store256d012(&ptr2[j+0],r2);
            store256d012(&ptr3[j+0],r3);
        }
    }
    {
        int i =  N-3;
        double* ptr0 = &this->v_[(i+0)*this->alignN];
        double* ptr1 = &this->v_[(i+1)*this->alignN];
        double* ptr2 = &this->v_[(i+2)*this->alignN];
        __m256d maxVal0 = _mm256_setzero_pd();
        __m256d maxVal1 = _mm256_setzero_pd();
        __m256d maxVal2 = _mm256_setzero_pd();
        const __m256d absmask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
        for(int j = 0; j < this->N-3; j=j+4)
        {
            __m256d absVal0 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
            __m256d absVal2 = _mm256_and_pd(load256d(&ptr2[j+0]), absmask);
            maxVal2 = _mm256_max_pd(absVal2,maxVal2);
        }
        {
            int j = this->N-3;
            __m256d absVal0 = _mm256_and_pd(load256d012(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d012(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
            __m256d absVal2 = _mm256_and_pd(load256d012(&ptr2[j+0]), absmask);
            maxVal2 = _mm256_max_pd(absVal2,maxVal2);
        }

        maxVal0 = broadcastMax(maxVal0);
        maxVal1 = broadcastMax(maxVal1);
        maxVal2 = broadcastMax(maxVal2);

        __m256d temp1 = _mm256_unpacklo_pd(maxVal0, maxVal1);
        __m256d temp2 = _mm256_unpacklo_pd(maxVal2, maxVal2);
        __m256d result = _mm256_permute2f128_pd(temp1, temp2, 0x20);
        __m256d ones = _mm256_set1_pd(1.0);
        result =  _mm256_div_pd(ones, result);
        storeu256d012(&this->rowScales[i+0],result);

        __m256d recipMaxVal0 = _mm256_permute4x64_pd(result, 0b00000000);
        __m256d recipMaxVal1 = _mm256_permute4x64_pd(result, 0b01010101);
        __m256d recipMaxVal2 = _mm256_permute4x64_pd(result, 0b10101010);
        for(int j = 0; j < this->N-3; j=j+4)
        {
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d(&ptr2[j+0]),recipMaxVal2);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
            store256d(&ptr2[j+0],r2);
        }
        {
            int j = this->N-3;
            __m256d r0 = mul256d(load256d012(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d012(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d012(&ptr2[j+0]),recipMaxVal2);
            store256d012(&ptr0[j+0],r0);
            store256d012(&ptr1[j+0],r1);
            store256d012(&ptr2[j+0],r2);
        }
    }
}

// Remain == 2 variant of the row equilibration.
void FastChemistry::LUsolver::scaleA2
(
)
{
    int remain8 = N%8;
    for(int i = 0; i < N-2; i=i+4)
    {
        double* ptr0 = &this->v_[(i+0)*this->alignN];
        double* ptr1 = &this->v_[(i+1)*this->alignN];
        double* ptr2 = &this->v_[(i+2)*this->alignN];
        double* ptr3 = &this->v_[(i+3)*this->alignN];

        __m256d maxVal0 = _mm256_setzero_pd();
        __m256d maxVal1 = _mm256_setzero_pd();
        __m256d maxVal2 = _mm256_setzero_pd();
        __m256d maxVal3 = _mm256_setzero_pd();
        const __m256d absmask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
        for(int j = 0; j < this->N-remain8; j=j+8)
        {

            __m256d absVal00 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            __m256d absVal04 = _mm256_and_pd(load256d(&ptr0[j+4]), absmask);
            maxVal0 = _mm256_max_pd(absVal00,maxVal0);
            maxVal0 = _mm256_max_pd(absVal04,maxVal0);
            __m256d absVal10 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            __m256d absVal14 = _mm256_and_pd(load256d(&ptr1[j+4]), absmask);
            maxVal1 = _mm256_max_pd(absVal10,maxVal1);
            maxVal1 = _mm256_max_pd(absVal14,maxVal1);
            __m256d absVal20 = _mm256_and_pd(load256d(&ptr2[j+0]), absmask);
            __m256d absVal24 = _mm256_and_pd(load256d(&ptr2[j+4]), absmask);
            maxVal2 = _mm256_max_pd(absVal20,maxVal2);
            maxVal2 = _mm256_max_pd(absVal24,maxVal2);
            __m256d absVal30 = _mm256_and_pd(load256d(&ptr3[j+0]), absmask);
            __m256d absVal34 = _mm256_and_pd(load256d(&ptr3[j+4]), absmask);
            maxVal3 = _mm256_max_pd(absVal30,maxVal3);
            maxVal3 = _mm256_max_pd(absVal34,maxVal3);
        }
        if(remain8==6)
        {
            int j = this->N-6;
            __m256d absVal0 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
            __m256d absVal2 = _mm256_and_pd(load256d(&ptr2[j+0]), absmask);
            maxVal2 = _mm256_max_pd(absVal2,maxVal2);
            __m256d absVal3 = _mm256_and_pd(load256d(&ptr3[j+0]), absmask);
            maxVal3 = _mm256_max_pd(absVal3,maxVal3);
        }
        {
            int j = this->N-2;
            __m256d absVal0 = _mm256_and_pd(load256d01(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d01(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
            __m256d absVal2 = _mm256_and_pd(load256d01(&ptr2[j+0]), absmask);
            maxVal2 = _mm256_max_pd(absVal2,maxVal2);
            __m256d absVal3 = _mm256_and_pd(load256d01(&ptr3[j+0]), absmask);
            maxVal3 = _mm256_max_pd(absVal3,maxVal3);
        }

        maxVal0 = broadcastMax(maxVal0);
        maxVal1 = broadcastMax(maxVal1);
        maxVal2 = broadcastMax(maxVal2);
        maxVal3 = broadcastMax(maxVal3);

        __m256d temp1 = _mm256_unpacklo_pd(maxVal0, maxVal1);
        __m256d temp2 = _mm256_unpacklo_pd(maxVal2, maxVal3);
        __m256d result = _mm256_permute2f128_pd(temp1, temp2, 0x20);
        __m256d ones = _mm256_set1_pd(1.0);
        result =  _mm256_div_pd(ones, result);
        storeu256d(&this->rowScales[i+0],result);
        __m256d recipMaxVal0 = _mm256_permute4x64_pd(result, 0b00000000);
        __m256d recipMaxVal1 = _mm256_permute4x64_pd(result, 0b01010101);
        __m256d recipMaxVal2 = _mm256_permute4x64_pd(result, 0b10101010);
        __m256d recipMaxVal3 = _mm256_permute4x64_pd(result, 0b11111111);

        for(int j = 0; j < this->N-remain8; j=j+8)
        {
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d(&ptr3[j+0]),recipMaxVal3);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
            store256d(&ptr2[j+0],r2);
            store256d(&ptr3[j+0],r3);
            __m256d r4 = mul256d(load256d(&ptr0[j+4]),recipMaxVal0);
            __m256d r5 = mul256d(load256d(&ptr1[j+4]),recipMaxVal1);
            __m256d r6 = mul256d(load256d(&ptr2[j+4]),recipMaxVal2);
            __m256d r7 = mul256d(load256d(&ptr3[j+4]),recipMaxVal3);
            store256d(&ptr0[j+4],r4);
            store256d(&ptr1[j+4],r5);
            store256d(&ptr2[j+4],r6);
            store256d(&ptr3[j+4],r7);
        }
        if(remain8==6)
        {
            int j = this->N-6;
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d(&ptr3[j+0]),recipMaxVal3);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
            store256d(&ptr2[j+0],r2);
            store256d(&ptr3[j+0],r3);
        }
        {

            int j = this->N-2;
            __m256d r0 = mul256d(load256d01(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d01(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d01(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d01(&ptr3[j+0]),recipMaxVal3);
            store256d01(&ptr0[j+0],r0);
            store256d01(&ptr1[j+0],r1);
            store256d01(&ptr2[j+0],r2);
            store256d01(&ptr3[j+0],r3);
        }
    }
    {
        int i =  N-2;
        double* ptr0 = &this->v_[(i+0)*this->alignN];
        double* ptr1 = &this->v_[(i+1)*this->alignN];
        __m256d maxVal0 = _mm256_setzero_pd();
        __m256d maxVal1 = _mm256_setzero_pd();
        const __m256d absmask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
        for(int j = 0; j < this->N-2; j=j+4)
        {
            __m256d absVal0 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
        }
        {
            int j = this->N-2;
            __m256d absVal0 = _mm256_and_pd(load256d01(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d01(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
        }

        maxVal0 = broadcastMax(maxVal0);
        maxVal1 = broadcastMax(maxVal1);


        __m256d temp1 = _mm256_unpacklo_pd(maxVal0, maxVal1);
        __m256d temp2 = _mm256_unpacklo_pd(maxVal0, maxVal1);
        __m256d result = _mm256_permute2f128_pd(temp1, temp2, 0x20);
        __m256d ones = _mm256_set1_pd(1.0);
        result =  _mm256_div_pd(ones, result);
        storeu256d01(&this->rowScales[i+0],result);

        __m256d recipMaxVal0 = _mm256_permute4x64_pd(result, 0b00000000);
        __m256d recipMaxVal1 = _mm256_permute4x64_pd(result, 0b01010101);

        for(int j = 0; j < this->N-2; j=j+4)
        {
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
        }
        {
            int j = this->N-2;
            __m256d r0 = mul256d(load256d01(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d01(&ptr1[j+0]),recipMaxVal1);
            store256d01(&ptr0[j+0],r0);
            store256d01(&ptr1[j+0],r1);
        }
    }
}

// Remain == 1 variant of the row equilibration.
void FastChemistry::LUsolver::scaleA1
(
)
{
    int remain8 = N%8;
    for(int i = 0; i < N-1; i=i+4)
    {
        double* ptr0 = &this->v_[(i+0)*this->alignN];
        double* ptr1 = &this->v_[(i+1)*this->alignN];
        double* ptr2 = &this->v_[(i+2)*this->alignN];
        double* ptr3 = &this->v_[(i+3)*this->alignN];

        __m256d maxVal0 = _mm256_setzero_pd();
        __m256d maxVal1 = _mm256_setzero_pd();
        __m256d maxVal2 = _mm256_setzero_pd();
        __m256d maxVal3 = _mm256_setzero_pd();
        const __m256d absmask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
        for(int j = 0; j < this->N-remain8; j=j+8)
        {

            __m256d absVal00 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            __m256d absVal04 = _mm256_and_pd(load256d(&ptr0[j+4]), absmask);
            maxVal0 = _mm256_max_pd(absVal00,maxVal0);
            maxVal0 = _mm256_max_pd(absVal04,maxVal0);
            __m256d absVal10 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            __m256d absVal14 = _mm256_and_pd(load256d(&ptr1[j+4]), absmask);
            maxVal1 = _mm256_max_pd(absVal10,maxVal1);
            maxVal1 = _mm256_max_pd(absVal14,maxVal1);
            __m256d absVal20 = _mm256_and_pd(load256d(&ptr2[j+0]), absmask);
            __m256d absVal24 = _mm256_and_pd(load256d(&ptr2[j+4]), absmask);
            maxVal2 = _mm256_max_pd(absVal20,maxVal2);
            maxVal2 = _mm256_max_pd(absVal24,maxVal2);
            __m256d absVal30 = _mm256_and_pd(load256d(&ptr3[j+0]), absmask);
            __m256d absVal34 = _mm256_and_pd(load256d(&ptr3[j+4]), absmask);
            maxVal3 = _mm256_max_pd(absVal30,maxVal3);
            maxVal3 = _mm256_max_pd(absVal34,maxVal3);
        }
        if(remain8==5)
        {
            int j = this->N-5;
            __m256d absVal0 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
            __m256d absVal2 = _mm256_and_pd(load256d(&ptr2[j+0]), absmask);
            maxVal2 = _mm256_max_pd(absVal2,maxVal2);
            __m256d absVal3 = _mm256_and_pd(load256d(&ptr3[j+0]), absmask);
            maxVal3 = _mm256_max_pd(absVal3,maxVal3);
        }
        {
            int j = this->N-1;
            __m256d absVal0 = _mm256_and_pd(load256d01(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d01(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
            __m256d absVal2 = _mm256_and_pd(load256d01(&ptr2[j+0]), absmask);
            maxVal2 = _mm256_max_pd(absVal2,maxVal2);
            __m256d absVal3 = _mm256_and_pd(load256d01(&ptr3[j+0]), absmask);
            maxVal3 = _mm256_max_pd(absVal3,maxVal3);
        }

        maxVal0 = broadcastMax(maxVal0);
        maxVal1 = broadcastMax(maxVal1);
        maxVal2 = broadcastMax(maxVal2);
        maxVal3 = broadcastMax(maxVal3);

        __m256d temp1 = _mm256_unpacklo_pd(maxVal0, maxVal1);
        __m256d temp2 = _mm256_unpacklo_pd(maxVal2, maxVal3);
        __m256d result = _mm256_permute2f128_pd(temp1, temp2, 0x20);
        __m256d ones = _mm256_set1_pd(1.0);
        result =  _mm256_div_pd(ones, result);
        storeu256d(&this->rowScales[i+0],result);
        __m256d recipMaxVal0 = _mm256_permute4x64_pd(result, 0b00000000);
        __m256d recipMaxVal1 = _mm256_permute4x64_pd(result, 0b01010101);
        __m256d recipMaxVal2 = _mm256_permute4x64_pd(result, 0b10101010);
        __m256d recipMaxVal3 = _mm256_permute4x64_pd(result, 0b11111111);

        for(int j = 0; j < this->N-remain8; j=j+8)
        {
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d(&ptr3[j+0]),recipMaxVal3);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
            store256d(&ptr2[j+0],r2);
            store256d(&ptr3[j+0],r3);
            __m256d r4 = mul256d(load256d(&ptr0[j+4]),recipMaxVal0);
            __m256d r5 = mul256d(load256d(&ptr1[j+4]),recipMaxVal1);
            __m256d r6 = mul256d(load256d(&ptr2[j+4]),recipMaxVal2);
            __m256d r7 = mul256d(load256d(&ptr3[j+4]),recipMaxVal3);
            store256d(&ptr0[j+4],r4);
            store256d(&ptr1[j+4],r5);
            store256d(&ptr2[j+4],r6);
            store256d(&ptr3[j+4],r7);
        }
        if(remain8==5)
        {
            int j = this->N-5;
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d(&ptr3[j+0]),recipMaxVal3);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
            store256d(&ptr2[j+0],r2);
            store256d(&ptr3[j+0],r3);
        }
        {

            int j = this->N-1;
            __m256d r0 = mul256d(load256d0(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d0(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d0(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d0(&ptr3[j+0]),recipMaxVal3);
            store256d0(&ptr0[j+0],r0);
            store256d0(&ptr1[j+0],r1);
            store256d0(&ptr2[j+0],r2);
            store256d0(&ptr3[j+0],r3);
        }
    }
    {
        int i =  N-1;
        double* ptr0 = &this->v_[(i+0)*this->alignN];
        __m256d maxVal0 = _mm256_setzero_pd();
        const __m256d absmask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
        for(int j = 0; j < this->N-1; j=j+4)
        {
            __m256d absVal0 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
        }
        {
            int j = this->N-1;
            __m256d absVal0 = _mm256_and_pd(load256d01(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
        }

        maxVal0 = broadcastMax(maxVal0);

        __m256d ones = _mm256_set1_pd(1.0);
        __m256d result =  _mm256_div_pd(ones, maxVal0);
        storeu256d0(&this->rowScales[i+0],result);

        __m256d recipMaxVal0 = _mm256_permute4x64_pd(result, 0b00000000);

        for(int j = 0; j < this->N-1; j=j+4)
        {
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            store256d(&ptr0[j+0],r0);
        }
        {
            int j = this->N-1;
            __m256d r0 = mul256d(load256d0(&ptr0[j+0]),recipMaxVal0);
            store256d0(&ptr0[j+0],r0);
        }
    }
}

// Remain == 0 variant of the row equilibration.
void FastChemistry::LUsolver::scaleA0
(
)
{
    int remain8 = N%8;
    for(int i = 0; i < N; i=i+4)
    {
        double* ptr0 = &this->v_[(i+0)*this->alignN];
        double* ptr1 = &this->v_[(i+1)*this->alignN];
        double* ptr2 = &this->v_[(i+2)*this->alignN];
        double* ptr3 = &this->v_[(i+3)*this->alignN];

        __m256d maxVal0 = _mm256_setzero_pd();
        __m256d maxVal1 = _mm256_setzero_pd();
        __m256d maxVal2 = _mm256_setzero_pd();
        __m256d maxVal3 = _mm256_setzero_pd();
        const __m256d absmask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
        for(int j = 0; j < this->N-remain8; j=j+8)
        {

            __m256d absVal00 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            __m256d absVal04 = _mm256_and_pd(load256d(&ptr0[j+4]), absmask);
            maxVal0 = _mm256_max_pd(absVal00,maxVal0);
            maxVal0 = _mm256_max_pd(absVal04,maxVal0);
            __m256d absVal10 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            __m256d absVal14 = _mm256_and_pd(load256d(&ptr1[j+4]), absmask);
            maxVal1 = _mm256_max_pd(absVal10,maxVal1);
            maxVal1 = _mm256_max_pd(absVal14,maxVal1);
            __m256d absVal20 = _mm256_and_pd(load256d(&ptr2[j+0]), absmask);
            __m256d absVal24 = _mm256_and_pd(load256d(&ptr2[j+4]), absmask);
            maxVal2 = _mm256_max_pd(absVal20,maxVal2);
            maxVal2 = _mm256_max_pd(absVal24,maxVal2);
            __m256d absVal30 = _mm256_and_pd(load256d(&ptr3[j+0]), absmask);
            __m256d absVal34 = _mm256_and_pd(load256d(&ptr3[j+4]), absmask);
            maxVal3 = _mm256_max_pd(absVal30,maxVal3);
            maxVal3 = _mm256_max_pd(absVal34,maxVal3);
        }
        if(remain8==4)
        {
            int j = this->N-4;
            __m256d absVal0 = _mm256_and_pd(load256d(&ptr0[j+0]), absmask);
            maxVal0 = _mm256_max_pd(absVal0,maxVal0);
            __m256d absVal1 = _mm256_and_pd(load256d(&ptr1[j+0]), absmask);
            maxVal1 = _mm256_max_pd(absVal1,maxVal1);
            __m256d absVal2 = _mm256_and_pd(load256d(&ptr2[j+0]), absmask);
            maxVal2 = _mm256_max_pd(absVal2,maxVal2);
            __m256d absVal3 = _mm256_and_pd(load256d(&ptr3[j+0]), absmask);
            maxVal3 = _mm256_max_pd(absVal3,maxVal3);
        }

        maxVal0 = broadcastMax(maxVal0);
        maxVal1 = broadcastMax(maxVal1);
        maxVal2 = broadcastMax(maxVal2);
        maxVal3 = broadcastMax(maxVal3);

        __m256d temp1 = _mm256_unpacklo_pd(maxVal0, maxVal1);
        __m256d temp2 = _mm256_unpacklo_pd(maxVal2, maxVal3);
        __m256d result = _mm256_permute2f128_pd(temp1, temp2, 0x20);
        __m256d ones = _mm256_set1_pd(1.0);
        result =  _mm256_div_pd(ones, result);
        storeu256d(&this->rowScales[i+0],result);
        __m256d recipMaxVal0 = _mm256_permute4x64_pd(result, 0b00000000);
        __m256d recipMaxVal1 = _mm256_permute4x64_pd(result, 0b01010101);
        __m256d recipMaxVal2 = _mm256_permute4x64_pd(result, 0b10101010);
        __m256d recipMaxVal3 = _mm256_permute4x64_pd(result, 0b11111111);

        for(int j = 0; j < this->N-remain8; j=j+8)
        {
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d(&ptr3[j+0]),recipMaxVal3);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
            store256d(&ptr2[j+0],r2);
            store256d(&ptr3[j+0],r3);
            __m256d r4 = mul256d(load256d(&ptr0[j+4]),recipMaxVal0);
            __m256d r5 = mul256d(load256d(&ptr1[j+4]),recipMaxVal1);
            __m256d r6 = mul256d(load256d(&ptr2[j+4]),recipMaxVal2);
            __m256d r7 = mul256d(load256d(&ptr3[j+4]),recipMaxVal3);
            store256d(&ptr0[j+4],r4);
            store256d(&ptr1[j+4],r5);
            store256d(&ptr2[j+4],r6);
            store256d(&ptr3[j+4],r7);
        }
        if(remain8==4)
        {
            int j = this->N-4;
            __m256d r0 = mul256d(load256d(&ptr0[j+0]),recipMaxVal0);
            __m256d r1 = mul256d(load256d(&ptr1[j+0]),recipMaxVal1);
            __m256d r2 = mul256d(load256d(&ptr2[j+0]),recipMaxVal2);
            __m256d r3 = mul256d(load256d(&ptr3[j+0]),recipMaxVal3);
            store256d(&ptr0[j+0],r0);
            store256d(&ptr1[j+0],r1);
            store256d(&ptr2[j+0],r2);
            store256d(&ptr3[j+0],r3);
        }
    }
}

// Apply the stored row scales to the right-hand side: b[i] *= rowScales[i].
void FastChemistry::LUsolver::scaleb
(
    double* __restrict__ b
)
{
    for(int i = 0; i < this->N; i++)
    {
        b[i] *= this->rowScales[i];
    }
}

// Full block 4x4 LU decomposition: row equilibration, block LU, pivoting,
// substitution and trailing updates (calls the Remain-specific kernels).
void FastChemistry::LUsolver::Block4LUDecompose
(
)
{

    int remain = this->N%4;
    for(int i = 0; i < this->N;i++)
    {
        this->pivotIndice_[i] = i;
    }
    int times = (this->N-remain)/4;
    if(remain==0)
    {
        this->scaleA0();
        for(int i = 0; i < times; i ++)
        {
            for(int j = 0; j <N;j++)
            {
                rowPtr[j] = &v_[j*alignN];
            }
            int k0 = i*4;
            int k1 = i*4+4;
            this->LUDecompose4(k0);
            this->forwardSubstitute4_0(k0,k1);
            this->permutation0(k0,k1);
            this->backSubstitute4_0(k0,k1);
            this->UpdateL22U22_Vec_38(k0,k1);
        }
    }
    else if(remain==1)
    {
        this->scaleA1();
        for(int i = 0; i < times; i ++)
        {
            for(int j = 0; j <N;j++)
            {
                rowPtr[j] = &v_[j*alignN];
            }
            int k0 = i*4;
            int k1 = i*4+4;
            this->LUDecompose4(k0);
            this->forwardSubstitute4_1(k0,k1);
            this->permutation1(k0,k1);
            this->backSubstitute4_1(k0,k1);
            this->UpdateL22U22_Vec_38(k0,k1);
        }
        if(v_[(N-1)*alignN+N-1]==0)
        {
            v_[(N-1)*alignN+N-1] = FastChemistry::LuLimiter;
        }
        this->invD[N-1] = 1.0/v_[(N-1)*alignN+N-1];

    }
    else if(remain==2)
    {
        this->scaleA2();
        for(int i = 0; i < times; i ++)
        {
            for(int j = 0; j <N;j++)
            {
                rowPtr[j] = &v_[j*alignN];
            }
            int k0 = i*4;
            int k1 = i*4+4;
            this->LUDecompose4(k0);
            this->forwardSubstitute4_2(k0,k1);
            this->permutation2(k0,k1);
            this->backSubstitute4_2(k0,k1);
            this->UpdateL22U22_Vec_38(k0,k1);
        }
        this->LUDecompose4_2();
        this->permutation_2();
    }
    else
    {
        this->scaleA3();
        for(int i = 0; i < times; i ++)
        {
            for(int j = 0; j <N;j++)
            {
                rowPtr[j] = &v_[j*alignN];
            }
            int k0 = i*4;
            int k1 = i*4+4;
            this->LUDecompose4(k0);
            this->forwardSubstitute4_3(k0,k1);
            this->permutation3(k0,k1);
            this->backSubstitute4_3(k0,k1);
            this->UpdateL22U22_Vec_38(k0,k1);
        }
        this->LUDecompose4_3();
        this->permutation_3();
    }
}

// Reassign the matrix data pointer.
void FastChemistry::LUsolver::ReAssign(double* externalData)
{this->v_ = externalData;}
