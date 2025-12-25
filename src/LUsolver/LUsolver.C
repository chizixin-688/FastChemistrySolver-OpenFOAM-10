/*---------------------------------------------------------------------------*\
  Description
      Implementation file for LU decomposition of matrix. 

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <iostream>
#include <cmath>

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
    this->owner = false;

    this->N = size;
    //this->alignN = size+(4-(size%4));
    this->alignN = ((size+3)/4)*4;
    this->v_ = externalData;
    this->Remain=this->N%4;


    this->pivotIndice_ = new unsigned int[N];
    for(unsigned int i = 0; i< N; i++)
    {
        this->pivotIndice_[i] = i;
    }
    rowPtr.resize(N);
    for(size_t i = 0; i <N;i++)
    {
        rowPtr[i] = &v_[i*alignN];
    }
    invD.resize(N);
}

FastChemistry::LUsolver::~LUsolver()
{
    delete[] this->pivotIndice_;
    this->pivotIndice_ = nullptr;
    if(this->owner==false)
    {
        this->v_ = nullptr;        
    }
    else
    {
        free(this->v_);
    }
}


void FastChemistry::LUsolver::printMatrix(double*A,int mRows, int nCol)
{
    for(int i = 0;i<mRows;i++)
    {
        for(int j =0;j<nCol-1;j++)
        {
            std::cout<<A[i*nCol+j]<<" ";
        }
        std::cout<<A[i*nCol+nCol-1]<<std::endl;
    }
}

void FastChemistry::LUsolver::printMatrix()
{
    for(unsigned int i = 0;i<this->N;i++)
    {
        for(unsigned int j =0;j<this->N-1;j++)
        {
            std::cout<<this->v_[i*this->N+j]<<" ";
        }
        std::cout<<this->v_[i*this->N+this->N-1]<<std::endl;
    }
}

void FastChemistry::LUsolver::printPivotIndice()
{
    for(unsigned int i = 0;i<this->N;i++)
    {
        std::cout<<this->pivotIndice_[i]<<std::endl;
    }    
}

void FastChemistry::LUsolver::LUDecompose4
(
    unsigned int k0
)
{

// step 1
    unsigned int iMax = 0;
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
        temp = rowk2[k0+1];
    }

    if (temp<std::fabs(rowk3[k0+1]))
    {
        iMax = 3;
        temp = rowk3[k0+1];        
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

void FastChemistry::LUsolver::LUDecompose4_3
(
)
{
    double Array0[4];
    {
        int k0 = this->N - 3;

// step 1
        unsigned int iMax = 0;
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

void FastChemistry::LUsolver::forwardSitituate4_0
(
    unsigned int k0,
    unsigned int k1
)
{
    double* __restrict__ rowK0 = rowPtr[k0+0];
    double* __restrict__ rowK1 = rowPtr[k0+1];
    double* __restrict__ rowK2 = rowPtr[k0+2];
    double* __restrict__ rowK3 = rowPtr[k0+3];

    const double L10 = rowK1[k0+0];
    const double L20 = rowK2[k0+0];
    const double L30 = rowK3[k0+0];
    const double L21 = rowK2[k0+1];
    const double L31 = rowK3[k0+1];
    const double L32 = rowK3[k0+2];
    __m256d negL10v = _mm256_set1_pd(-L10);
    __m256d negL20v = _mm256_set1_pd(-L20);
    __m256d negL21v = _mm256_set1_pd(-L21);
    __m256d negL30v = _mm256_set1_pd(-L30);
    __m256d negL31v = _mm256_set1_pd(-L31);
    __m256d negL32v = _mm256_set1_pd(-L32);

    for(unsigned int i = k1; i <= this->N-4; i=i+4)
    {
        __m256d U00v = load256d(&rowK0[i]);   

        __m256d U10v = load256d(&rowK1[i]);
        U10v = _mm256_fmadd_pd(negL10v,U00v,U10v);
        store256d(&rowK1[i],U10v);

        __m256d U20v = load256d(&rowK2[i]);
        U20v = _mm256_fmadd_pd(negL20v,U00v,U20v);
        U20v = _mm256_fmadd_pd(negL21v,U10v,U20v);
        store256d(&rowK2[i],U20v);            
            
        __m256d U30v = load256d(&rowK3[i]);
        U30v = _mm256_fmadd_pd(negL30v,U00v,U30v);
        U30v = _mm256_fmadd_pd(negL31v,U10v,U30v);
        U30v = _mm256_fmadd_pd(negL32v,U20v,U30v);
        store256d(&rowK3[i],U30v);  
    }
}

void FastChemistry::LUsolver::forwardSitituate4_1
(
    unsigned int k0,
    unsigned int k1
)
{
    double* __restrict__ rowK0 = rowPtr[k0+0];
    double* __restrict__ rowK1 = rowPtr[k0+1];
    double* __restrict__ rowK2 = rowPtr[k0+2];
    double* __restrict__ rowK3 = rowPtr[k0+3];

    const double L10 = rowK1[k0+0];
    const double L20 = rowK2[k0+0];
    const double L30 = rowK3[k0+0];
    const double L21 = rowK2[k0+1];
    const double L31 = rowK3[k0+1];
    const double L32 = rowK3[k0+2];
    __m256d negL10v = _mm256_set1_pd(-L10);
    __m256d negL20v = _mm256_set1_pd(-L20);
    __m256d negL21v = _mm256_set1_pd(-L21);
    __m256d negL30v = _mm256_set1_pd(-L30);
    __m256d negL31v = _mm256_set1_pd(-L31);
    __m256d negL32v = _mm256_set1_pd(-L32);

    for(unsigned int i = k1; i <= this->N-4; i=i+4)
    {
        __m256d U00v = load256d(&rowK0[i]);   

        __m256d U10v = load256d(&rowK1[i]);
        U10v = _mm256_fmadd_pd(negL10v,U00v,U10v);
        store256d(&rowK1[i],U10v);

        __m256d U20v = load256d(&rowK2[i]);
        U20v = _mm256_fmadd_pd(negL20v,U00v,U20v);
        U20v = _mm256_fmadd_pd(negL21v,U10v,U20v);
        store256d(&rowK2[i],U20v);            

        __m256d U30v = load256d(&rowK3[i]);
        U30v = _mm256_fmadd_pd(negL30v,U00v,U30v);
        U30v = _mm256_fmadd_pd(negL31v,U10v,U30v);
        U30v = _mm256_fmadd_pd(negL32v,U20v,U30v);
        store256d(&rowK3[i],U30v);  
    }
    {
        {
            rowK1[N-1] = rowK1[N-1] 
                         - L10*rowK0[N-1];
            rowK2[N-1] = rowK2[N-1] 
                         - L20*rowK0[N-1]
                         - L21*rowK1[N-1];         
            rowK3[N-1] = rowK3[N-1] 
                         - L30*rowK0[N-1]
                         - L31*rowK1[N-1]
                         - L32*rowK2[N-1];
        }
               
    }
}

void FastChemistry::LUsolver::forwardSitituate4_2
(
    unsigned int k0,
    unsigned int k1
)
{
    double* __restrict__ rowK0 = rowPtr[k0+0];
    double* __restrict__ rowK1 = rowPtr[k0+1];
    double* __restrict__ rowK2 = rowPtr[k0+2];
    double* __restrict__ rowK3 = rowPtr[k0+3];

    const double L10 = rowK1[k0+0];
    const double L20 = rowK2[k0+0];
    const double L30 = rowK3[k0+0];
    const double L21 = rowK2[k0+1];
    const double L31 = rowK3[k0+1];
    const double L32 = rowK3[k0+2];
    __m256d negL10v = _mm256_set1_pd(-L10);
    __m256d negL20v = _mm256_set1_pd(-L20);
    __m256d negL21v = _mm256_set1_pd(-L21);
    __m256d negL30v = _mm256_set1_pd(-L30);
    __m256d negL31v = _mm256_set1_pd(-L31);
    __m256d negL32v = _mm256_set1_pd(-L32);

    for(unsigned int i = k1; i <= this->N-4; i=i+4)
    {
        __m256d U00v = load256d(&rowK0[i]);   

        __m256d U10v = load256d(&rowK1[i]);
        U10v = _mm256_fmadd_pd(negL10v,U00v,U10v);
        store256d(&rowK1[i],U10v);

        __m256d U20v = load256d(&rowK2[i]);
        U20v = _mm256_fmadd_pd(negL20v,U00v,U20v);
        U20v = _mm256_fmadd_pd(negL21v,U10v,U20v);
        store256d(&rowK2[i],U20v);            
            
        __m256d U30v = load256d(&rowK3[i]);
        U30v = _mm256_fmadd_pd(negL30v,U00v,U30v);
        U30v = _mm256_fmadd_pd(negL31v,U10v,U30v);
        U30v = _mm256_fmadd_pd(negL32v,U20v,U30v);
        store256d(&rowK3[i],U30v);  
    }

    {
        unsigned int i = N-2;
            rowK1[i+0] = rowK1[i+0] //U10 = A10 -L10*U00
                         - L10*rowK0[i+0];
            rowK1[i+1] = rowK1[i+1] //U11 = A11 -L10*U01
                         - L10*rowK0[i+1];
            rowK2[i+0] = rowK2[i+0] //U20 = A20 - L20*U00 - L21*U10
                         - L20*rowK0[i+0]
                         - L21*rowK1[i+0];
            rowK2[i+1] = rowK2[i+1] //U21 = A21 - L20*U01 - L21*U11
                         - L20*rowK0[i+1]
                         - L21*rowK1[i+1];         
            rowK3[i+0] = rowK3[i+0] //U30 = U30 - L30*U00 - L31*U10 - L32*U20
                         - L30*rowK0[i+0]
                         - L31*rowK1[i+0]
                         - L32*rowK2[i+0];
            rowK3[i+1] = rowK3[i+1] //U31 = U31 - L30*U01 - L31*U11 - L32*U21
                         - L30*rowK0[i+1]
                         - L31*rowK1[i+1]
                         - L32*rowK2[i+1];
        
    }

}

void FastChemistry::LUsolver::forwardSitituate4_3
(
    unsigned int k0,
    unsigned int k1
)
{
    double* __restrict__ rowK0 = rowPtr[k0+0];
    double* __restrict__ rowK1 = rowPtr[k0+1];
    double* __restrict__ rowK2 = rowPtr[k0+2];
    double* __restrict__ rowK3 = rowPtr[k0+3];

    const double L10 = rowK1[k0+0];
    const double L20 = rowK2[k0+0];
    const double L30 = rowK3[k0+0];
    const double L21 = rowK2[k0+1];
    const double L31 = rowK3[k0+1];
    const double L32 = rowK3[k0+2];
    __m256d negL10v = _mm256_set1_pd(-L10);
    __m256d negL20v = _mm256_set1_pd(-L20);
    __m256d negL21v = _mm256_set1_pd(-L21);
    __m256d negL30v = _mm256_set1_pd(-L30);
    __m256d negL31v = _mm256_set1_pd(-L31);
    __m256d negL32v = _mm256_set1_pd(-L32);

    for(unsigned int i = k1; i <= this->N-4; i=i+4)
    {
        __m256d U00v = load256d(&rowK0[i]);   

        __m256d U10v = load256d(&rowK1[i]);
        U10v = _mm256_fmadd_pd(negL10v,U00v,U10v);
        store256d(&rowK1[i],U10v);

        __m256d U20v = load256d(&rowK2[i]);
        U20v = _mm256_fmadd_pd(negL20v,U00v,U20v);
        U20v = _mm256_fmadd_pd(negL21v,U10v,U20v);
        store256d(&rowK2[i],U20v);            
            
        __m256d U30v = load256d(&rowK3[i]);
        U30v = _mm256_fmadd_pd(negL30v,U00v,U30v);
        U30v = _mm256_fmadd_pd(negL31v,U10v,U30v);
        U30v = _mm256_fmadd_pd(negL32v,U20v,U30v);
        store256d(&rowK3[i],U30v);  
    }

    {
        unsigned int i = N-3;
        {

            rowK1[i+0] = rowK1[i+0] //U10 = A10 -L10*U00
                         - L10*rowK0[i+0];
            rowK1[i+1] = rowK1[i+1] //U11 = A11 -L10*U01
                         - L10*rowK0[i+1];
            rowK1[i+2] = rowK1[i+2] //U12 = A12 -L10*U02
                         - L10*rowK0[i+2];

            rowK2[i+0] = rowK2[i+0] //U20 = A20 - L20*U00 - L21*U10
                         - L20*rowK0[i+0]
                         - L21*rowK1[i+0];
            rowK2[i+1] = rowK2[i+1] //U21 = A21 - L20*U01 - L21*U11
                         - L20*rowK0[i+1]
                         - L21*rowK1[i+1];
            rowK2[i+2] = rowK2[i+2] //U22 = A22 - L20*U02 - L21*U12
                         - L20*rowK0[i+2]
                         - L21*rowK1[i+2];         

            rowK3[i+0] = rowK3[i+0] //U30 = U30 - L30*U00 - L31*U10 - L32*U20
                         - L30*rowK0[i+0]
                         - L31*rowK1[i+0]
                         - L32*rowK2[i+0];
            rowK3[i+1] = rowK3[i+1] //U31 = U31 - L30*U01 - L31*U11 - L32*U21
                         - L30*rowK0[i+1]
                         - L31*rowK1[i+1]
                         - L32*rowK2[i+1];
            rowK3[i+2] = rowK3[i+2] //U32 = U32 - L30*U02 - L31*U12 - L32*U22
                         - L30*rowK0[i+2]
                         - L31*rowK1[i+2]
                         - L32*rowK2[i+2];
        } 
    }
}



void FastChemistry::LUsolver::backSitituate4_0
(
    unsigned int k0,
    unsigned int k1
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
    __m256d negU01v = _mm256_set1_pd(-U01);
    __m256d negU02v = _mm256_set1_pd(-U02);
    __m256d negU03v = _mm256_set1_pd(-U03);
    __m256d negU12v = _mm256_set1_pd(-U12);
    __m256d negU13v = _mm256_set1_pd(-U13);
    __m256d negU23v = _mm256_set1_pd(-U23);


    for(unsigned int i = k1; i <= this->N-4; i = i + 4)
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

        L01v = _mm256_fmadd_pd(L00v,negU01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = _mm256_fmadd_pd(L00v,negU02v,L02v);
        L02v = _mm256_fmadd_pd(L01v,negU12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = _mm256_fmadd_pd(L00v,negU03v,L03v);
        L03v = _mm256_fmadd_pd(L01v,negU13v,L03v);
        L03v = _mm256_fmadd_pd(L02v,negU23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);      

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);     
    }
}


void FastChemistry::LUsolver::backSitituate4_1
(
    unsigned int k0,
    unsigned int k1
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
    __m256d negU01v = _mm256_set1_pd(-U01);
    __m256d negU02v = _mm256_set1_pd(-U02);
    __m256d negU03v = _mm256_set1_pd(-U03);
    __m256d negU12v = _mm256_set1_pd(-U12);
    __m256d negU13v = _mm256_set1_pd(-U13);
    __m256d negU23v = _mm256_set1_pd(-U23);


    for(unsigned int i = k1; i <= this->N-4; i = i + 4)
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

        L01v = _mm256_fmadd_pd(L00v,negU01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = _mm256_fmadd_pd(L00v,negU02v,L02v);
        L02v = _mm256_fmadd_pd(L01v,negU12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = _mm256_fmadd_pd(L00v,negU03v,L03v);
        L03v = _mm256_fmadd_pd(L01v,negU13v,L03v);
        L03v = _mm256_fmadd_pd(L02v,negU23v,L03v);
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

void FastChemistry::LUsolver::backSitituate4_2
(
    unsigned int k0,
    unsigned int k1
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
    __m256d negU01v = _mm256_set1_pd(-U01);
    __m256d negU02v = _mm256_set1_pd(-U02);
    __m256d negU03v = _mm256_set1_pd(-U03);
    __m256d negU12v = _mm256_set1_pd(-U12);
    __m256d negU13v = _mm256_set1_pd(-U13);
    __m256d negU23v = _mm256_set1_pd(-U23);

    for(unsigned int i = k1; i <= this->N-4; i = i + 4)
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

        L01v = _mm256_fmadd_pd(L00v,negU01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = _mm256_fmadd_pd(L00v,negU02v,L02v);
        L02v = _mm256_fmadd_pd(L01v,negU12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = _mm256_fmadd_pd(L00v,negU03v,L03v);
        L03v = _mm256_fmadd_pd(L01v,negU13v,L03v);
        L03v = _mm256_fmadd_pd(L02v,negU23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);      

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);
    }

    {
        double* __restrict__ rowN2 = &v_[(N-2)*alignN];
        double* __restrict__ rowN1 = rowN2 + alignN;
        double& La0 = rowN2[k0+0];
        double& Lb0 = rowN1[k0+0];
        La0 = La0*invU00;
        Lb0 = Lb0*invU00;

        double& La1 = rowN2[k0+1];
        double& Lb1 = rowN1[k0+1];
        La1 = (La1 - La0*U01)*invU11;
        Lb1 = (Lb1 - Lb0*U01)*invU11;

        double& La2 = rowN2[k0+2];
        double& Lb2 = rowN1[k0+2];
        La2 = (La2 - La0*U02 - La1*U12)*invU22;
        Lb2 = (Lb2 - Lb0*U02 - Lb1*U12)*invU22;

        double& La3 = rowN2[k0+3];
        double& Lb3 = rowN1[k0+3];
        La3 = (La3 - La0*U03 - La1*U13 - La2*U23)*invU33;
        Lb3 = (Lb3 - Lb0*U03 - Lb1*U13 - Lb2*U23)*invU33;
    }
}

void FastChemistry::LUsolver::backSitituate4_3
(
    unsigned int k0,
    unsigned int k1
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
    __m256d negU01v = _mm256_set1_pd(-U01);
    __m256d negU02v = _mm256_set1_pd(-U02);
    __m256d negU03v = _mm256_set1_pd(-U03);
    __m256d negU12v = _mm256_set1_pd(-U12);
    __m256d negU13v = _mm256_set1_pd(-U13);
    __m256d negU23v = _mm256_set1_pd(-U23);

    for(unsigned int i = k1; i <= this->N-4; i = i + 4)
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

        L01v = _mm256_fmadd_pd(L00v,negU01v,L01v);
        L01v = _mm256_mul_pd(L01v,invU11v);

        L02v = _mm256_fmadd_pd(L00v,negU02v,L02v);
        L02v = _mm256_fmadd_pd(L01v,negU12v,L02v);
        L02v = _mm256_mul_pd(L02v,invU22v);

        L03v = _mm256_fmadd_pd(L00v,negU03v,L03v);
        L03v = _mm256_fmadd_pd(L01v,negU13v,L03v);
        L03v = _mm256_fmadd_pd(L02v,negU23v,L03v);
        L03v = _mm256_mul_pd(L03v,invU33v);      

        transpose4x4_pd(L00v,L01v,L02v,L03v);
        store256d(&rowi0[k0+0],L00v);
        store256d(&rowi1[k0+0],L01v);
        store256d(&rowi2[k0+0],L02v);
        store256d(&rowi3[k0+0],L03v);      
    }

    {
        double* __restrict__ rowN3 = &v_[(N-3)*alignN];
        double* __restrict__ rowN2 = rowN3+alignN;
        double* __restrict__ rowN1 = rowN2+alignN;

        double& La0 = rowN3[k0+0];
        double& Lb0 = rowN2[k0+0];
        double& Lc0 = rowN1[k0+0];
        La0 = La0*invU00;
        Lb0 = Lb0*invU00;
        Lc0 = Lc0*invU00;

        double& La1 = rowN3[k0+1];
        double& Lb1 = rowN2[k0+1];
        double& Lc1 = rowN1[k0+1];
        La1 = (La1 - La0*U01)*invU11;
        Lb1 = (Lb1 - Lb0*U01)*invU11;
        Lc1 = (Lc1 - Lc0*U01)*invU11;

        double& La2 = rowN3[k0+2];
        double& Lb2 = rowN2[k0+2];
        double& Lc2 = rowN1[k0+2];
        La2 = (La2 - La0*U02 - La1*U12)*invU22;
        Lb2 = (Lb2 - Lb0*U02 - Lb1*U12)*invU22;
        Lc2 = (Lc2 - Lc0*U02 - Lc1*U12)*invU22;

        double& La3 = rowN3[k0+3];
        double& Lb3 = rowN2[k0+3];
        double& Lc3 = rowN1[k0+3];
        La3 = (La3 - La0*U03 - La1*U13 - La2*U23)*invU33;
        Lb3 = (Lb3 - Lb0*U03 - Lb1*U13 - Lb2*U23)*invU33;
        Lc3 = (Lc3 - Lc0*U03 - Lc1*U13 - Lc2*U23)*invU33;
    }

}

void FastChemistry::LUsolver::UpdateL22U22_Vec2_0
(
    unsigned int k0,
    unsigned int k1
)
{
    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];    
    for(unsigned int i = k1; i < this->N; i=i+4)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];
        double* __restrict__ rowi3  = &v_[(i+3)*alignN];

        const double Li30 = rowi3[k0+0];
        const double Li31 = rowi3[k0+1];
        const double Li32 = rowi3[k0+2];
        const double Li33 = rowi3[k0+3];

        const double Li20 = rowi2[k0+0];
        const double Li21 = rowi2[k0+1];
        const double Li22 = rowi2[k0+2];
        const double Li23 = rowi2[k0+3];

        const double Li10 = rowi1[k0+0];
        const double Li11 = rowi1[k0+1];
        const double Li12 = rowi1[k0+2];
        const double Li13 = rowi1[k0+3];

        const double Li00 = rowi0[k0+0];
        const double Li01 = rowi0[k0+1];
        const double Li02 = rowi0[k0+2];
        const double Li03 = rowi0[k0+3];
        __m256d L00 =  _mm256_set1_pd( -Li00);
        __m256d L01 =  _mm256_set1_pd( -Li01);
        __m256d L02 =  _mm256_set1_pd( -Li02);
        __m256d L03 =  _mm256_set1_pd( -Li03);
        __m256d L10 =  _mm256_set1_pd( -Li10);
        __m256d L11 =  _mm256_set1_pd( -Li11);
        __m256d L12 =  _mm256_set1_pd( -Li12);
        __m256d L13 =  _mm256_set1_pd( -Li13);
        __m256d L20 =  _mm256_set1_pd( -Li20);
        __m256d L21 =  _mm256_set1_pd( -Li21);
        __m256d L22 =  _mm256_set1_pd( -Li22);
        __m256d L23 =  _mm256_set1_pd( -Li23);  
        __m256d L30 =  _mm256_set1_pd( -Li30);
        __m256d L31 =  _mm256_set1_pd( -Li31);
        __m256d L32 =  _mm256_set1_pd( -Li32);
        __m256d L33 =  _mm256_set1_pd( -Li33);            
        for(unsigned int j = k1; j < this->N; j=j+4)
        {
            __m256d U0x =  _mm256_loadu_pd(&rowk0[j]);
            __m256d U1x =  _mm256_loadu_pd(&rowk1[j]);
            __m256d U2x =  _mm256_loadu_pd(&rowk2[j]);
            __m256d U3x =  _mm256_loadu_pd(&rowk3[j]);

            __m256d A  =  _mm256_loadu_pd( &rowi0[j] );
            A = _mm256_fmadd_pd(L00,U0x,A);
            A = _mm256_fmadd_pd(L01,U1x,A);
            A = _mm256_fmadd_pd(L02,U2x,A);
            A = _mm256_fmadd_pd(L03,U3x,A);  
            _mm256_storeu_pd(&rowi0[j], A);

            A =  _mm256_loadu_pd( &rowi1[j] );
            A = _mm256_fmadd_pd(L10,U0x,A);
            A = _mm256_fmadd_pd(L11,U1x,A);
            A = _mm256_fmadd_pd(L12,U2x,A);
            A = _mm256_fmadd_pd(L13,U3x,A); 
            _mm256_storeu_pd(&rowi1[j], A);

            A =  _mm256_loadu_pd( &rowi2[j] );
            A = _mm256_fmadd_pd(L20,U0x,A);
            A = _mm256_fmadd_pd(L21,U1x,A);
            A = _mm256_fmadd_pd(L22,U2x,A);
            A = _mm256_fmadd_pd(L23,U3x,A); 
            _mm256_storeu_pd(&rowi2[j], A);

            A =  _mm256_loadu_pd( &rowi3[j] );
            A = _mm256_fmadd_pd(L30,U0x,A);
            A = _mm256_fmadd_pd(L31,U1x,A);
            A = _mm256_fmadd_pd(L32,U2x,A);
            A = _mm256_fmadd_pd(L33,U3x,A); 
            _mm256_storeu_pd(&rowi3[j], A);
        }
    }
}
void FastChemistry::LUsolver::UpdateL22U22_Vec2_1
(
    unsigned int k0,
    unsigned int k1
)
{
    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];        
    for(unsigned int i = k1; i < this->N -1; i=i+4)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];
        double* __restrict__ rowi3  = &v_[(i+3)*alignN]; 
        const double Li30 = rowi3[k0+0];
        const double Li31 = rowi3[k0+1];
        const double Li32 = rowi3[k0+2];
        const double Li33 = rowi3[k0+3];

        const double Li20 = rowi2[k0+0];
        const double Li21 = rowi2[k0+1];
        const double Li22 = rowi2[k0+2];
        const double Li23 = rowi2[k0+3];

        const double Li10 = rowi1[k0+0];
        const double Li11 = rowi1[k0+1];
        const double Li12 = rowi1[k0+2];
        const double Li13 = rowi1[k0+3];

        const double Li00 = rowi0[k0+0];
        const double Li01 = rowi0[k0+1];
        const double Li02 = rowi0[k0+2];
        const double Li03 = rowi0[k0+3];
        __m256d L00 =  _mm256_set1_pd(-Li00);
        __m256d L01 =  _mm256_set1_pd(-Li01);
        __m256d L02 =  _mm256_set1_pd(-Li02);
        __m256d L03 =  _mm256_set1_pd(-Li03);
        __m256d L10 =  _mm256_set1_pd(-Li10);
        __m256d L11 =  _mm256_set1_pd(-Li11);
        __m256d L12 =  _mm256_set1_pd(-Li12);
        __m256d L13 =  _mm256_set1_pd(-Li13);
        __m256d L20 =  _mm256_set1_pd(-Li20);
        __m256d L21 =  _mm256_set1_pd(-Li21);
        __m256d L22 =  _mm256_set1_pd(-Li22);
        __m256d L23 =  _mm256_set1_pd(-Li23);  
        __m256d L30 =  _mm256_set1_pd(-Li30);
        __m256d L31 =  _mm256_set1_pd(-Li31);
        __m256d L32 =  _mm256_set1_pd(-Li32);
        __m256d L33 =  _mm256_set1_pd(-Li33);                       
        for(unsigned int j = k1; j < this->N -1; j=j+4)
        {
            __m256d U0x =  _mm256_loadu_pd(&rowk0[j]);
            __m256d U1x =  _mm256_loadu_pd(&rowk1[j]);
            __m256d U2x =  _mm256_loadu_pd(&rowk2[j]);
            __m256d U3x =  _mm256_loadu_pd(&rowk3[j]);

            __m256d A  =  _mm256_loadu_pd( &rowi0[j] );
            A = _mm256_fmadd_pd(L00,U0x,A);
            A = _mm256_fmadd_pd(L01,U1x,A);
            A = _mm256_fmadd_pd(L02,U2x,A);
            A = _mm256_fmadd_pd(L03,U3x,A);  
            _mm256_storeu_pd(&rowi0[j], A);

            A =  _mm256_loadu_pd( &rowi1[j] );
            A = _mm256_fmadd_pd(L10,U0x,A);
            A = _mm256_fmadd_pd(L11,U1x,A);
            A = _mm256_fmadd_pd(L12,U2x,A);
            A = _mm256_fmadd_pd(L13,U3x,A); 
            _mm256_storeu_pd(&rowi1[j], A);

            A =  _mm256_loadu_pd( &rowi2[j] );
            A = _mm256_fmadd_pd(L20,U0x,A);
            A = _mm256_fmadd_pd(L21,U1x,A);
            A = _mm256_fmadd_pd(L22,U2x,A);
            A = _mm256_fmadd_pd(L23,U3x,A); 
            _mm256_storeu_pd(&rowi2[j], A);

            A =  _mm256_loadu_pd( &rowi3[j] );
            A = _mm256_fmadd_pd(L30,U0x,A);
            A = _mm256_fmadd_pd(L31,U1x,A);
            A = _mm256_fmadd_pd(L32,U2x,A);
            A = _mm256_fmadd_pd(L33,U3x,A); 
            _mm256_storeu_pd(&rowi3[j], A);                                 
        }
    }

    {
        double U0 = rowk0[N-1];
        double U1 = rowk1[N-1];
        double U2 = rowk2[N-1];
        double U3 = rowk3[N-1];  
        __m256d Uv = _mm256_setr_pd(U0,U1,U2,U3);

        for(unsigned int i = k1; i < this->N-1; i=i+4)
        {
            double* __restrict__ rowi0 = &v_[(i+0)*alignN];
            double* __restrict__ rowi1 = &v_[(i+1)*alignN];
            double* __restrict__ rowi2 = &v_[(i+2)*alignN];
            double* __restrict__ rowi3 = &v_[(i+3)*alignN];

            __m256d L0v = _mm256_loadu_pd(&rowi0[k0+0]);
            __m256d L1v = _mm256_loadu_pd(&rowi1[k0+0]);
            __m256d L2v = _mm256_loadu_pd(&rowi2[k0+0]);
            __m256d L3v = _mm256_loadu_pd(&rowi3[k0+0]); 

            __m256d r0 = _mm256_mul_pd(L0v,Uv);
            __m256d r1 = _mm256_mul_pd(L1v,Uv);
            __m256d r2 = _mm256_mul_pd(L2v,Uv);
            __m256d r3 = _mm256_mul_pd(L3v,Uv);

            __m256d r4 = hsum4x4(r0,r1,r2,r3);

            /*rowi0[N-1] = rowi0[N-1] - hsum4(_mm256_mul_pd(L0v,Uv));
            rowi1[N-1] = rowi1[N-1] - hsum4(_mm256_mul_pd(L1v,Uv));
            rowi2[N-1] = rowi2[N-1] - hsum4(_mm256_mul_pd(L2v,Uv));
            rowi3[N-1] = rowi3[N-1] - hsum4(_mm256_mul_pd(L3v,Uv));*/
            rowi0[N-1] = rowi0[N-1] - get_elem0(r4);
            rowi1[N-1] = rowi1[N-1] - get_elem1(r4);
            rowi2[N-1] = rowi2[N-1] - get_elem2(r4);
            rowi3[N-1] = rowi3[N-1] - get_elem3(r4);
        }

        double* __restrict__ rowN1 = &v_[(N-1)*alignN];
        double L0 = rowN1[k0+0];
        double L1 = rowN1[k0+1];
        double L2 = rowN1[k0+2];
        double L3 = rowN1[k0+3];
        __m256d L0v = _mm256_set1_pd(-L0);
        __m256d L1v = _mm256_set1_pd(-L1);
        __m256d L2v = _mm256_set1_pd(-L2);
        __m256d L3v = _mm256_set1_pd(-L3);
        __m256d Lv = _mm256_loadu_pd(&rowN1[k0+0]);
        for(unsigned int i = k1; i < this->N-1; i=i+4)  
        {
            __m256d rowN1v = _mm256_loadu_pd(&rowN1[i+0]);
            __m256d rowk0v = _mm256_loadu_pd(&rowk0[i+0]);
            __m256d rowk1v = _mm256_loadu_pd(&rowk1[i+0]);
            __m256d rowk2v = _mm256_loadu_pd(&rowk2[i+0]);
            __m256d rowk3v = _mm256_loadu_pd(&rowk3[i+0]);
            rowN1v = _mm256_fmadd_pd(rowk0v,L0v,rowN1v);
            rowN1v = _mm256_fmadd_pd(rowk1v,L1v,rowN1v);
            rowN1v = _mm256_fmadd_pd(rowk2v,L2v,rowN1v);
            rowN1v = _mm256_fmadd_pd(rowk3v,L3v,rowN1v);
            _mm256_storeu_pd(&rowN1[i+0],rowN1v);
        }
        {
            unsigned int i=N-1;
            double U0_ = rowk0[i];
            double U1_ = rowk1[i];
            double U2_ = rowk2[i];
            double U3_ = rowk3[i];
            __m256d temp = _mm256_setr_pd(U0_,U1_,U2_,U3_);
            rowN1[i] = rowN1[i] - hsum4(_mm256_mul_pd(Lv,temp));       
        }
    }

}


void FastChemistry::LUsolver::UpdateL22U22_Vec2_2
(
    unsigned int k0,
    unsigned int k1
)
{
    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];

    for(unsigned int i = k1; i < this->N -2; i=i+4)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];
        double* __restrict__ rowi3  = &v_[(i+3)*alignN];

        const double Li30 = rowi3[k0+0];
        const double Li31 = rowi3[k0+1];
        const double Li32 = rowi3[k0+2];
        const double Li33 = rowi3[k0+3];

        const double Li20 = rowi2[k0+0];
        const double Li21 = rowi2[k0+1];
        const double Li22 = rowi2[k0+2];
        const double Li23 = rowi2[k0+3];

        const double Li10 = rowi1[k0+0];
        const double Li11 = rowi1[k0+1];
        const double Li12 = rowi1[k0+2];
        const double Li13 = rowi1[k0+3];

        const double Li00 = rowi0[k0+0];
        const double Li01 = rowi0[k0+1];
        const double Li02 = rowi0[k0+2];
        const double Li03 = rowi0[k0+3];

        __m256d L00 =  _mm256_set1_pd(-Li00);
        __m256d L01 =  _mm256_set1_pd(-Li01);
        __m256d L02 =  _mm256_set1_pd(-Li02);
        __m256d L03 =  _mm256_set1_pd(-Li03);
        __m256d L10 =  _mm256_set1_pd(-Li10);
        __m256d L11 =  _mm256_set1_pd(-Li11);
        __m256d L12 =  _mm256_set1_pd(-Li12);
        __m256d L13 =  _mm256_set1_pd(-Li13);
        __m256d L20 =  _mm256_set1_pd(-Li20);
        __m256d L21 =  _mm256_set1_pd(-Li21);
        __m256d L22 =  _mm256_set1_pd(-Li22);
        __m256d L23 =  _mm256_set1_pd(-Li23);
        __m256d L30 =  _mm256_set1_pd(-Li30);
        __m256d L31 =  _mm256_set1_pd(-Li31);
        __m256d L32 =  _mm256_set1_pd(-Li32);
        __m256d L33 =  _mm256_set1_pd(-Li33);




        for(unsigned int j = k1; j < this->N -2; j=j+4)
        {
            __m256d U0x = _mm256_loadu_pd(&rowk0[j]);
            __m256d U1x = _mm256_loadu_pd(&rowk1[j]);
            __m256d U2x = _mm256_loadu_pd(&rowk2[j]);
            __m256d U3x = _mm256_loadu_pd(&rowk3[j]);

            __m256d A = _mm256_loadu_pd( &rowi0[j] );
            A = _mm256_fmadd_pd(L00,U0x,A);
            A = _mm256_fmadd_pd(L01,U1x,A);
            A = _mm256_fmadd_pd(L02,U2x,A);
            A = _mm256_fmadd_pd(L03,U3x,A);  
            _mm256_storeu_pd(&rowi0[j], A);

            A =  _mm256_loadu_pd( &rowi1[j] );
            A = _mm256_fmadd_pd(L10,U0x,A);
            A = _mm256_fmadd_pd(L11,U1x,A);
            A = _mm256_fmadd_pd(L12,U2x,A);
            A = _mm256_fmadd_pd(L13,U3x,A); 
            _mm256_storeu_pd(&rowi1[j], A);

            A =  _mm256_loadu_pd( &rowi2[j] );
            A = _mm256_fmadd_pd(L20,U0x,A);
            A = _mm256_fmadd_pd(L21,U1x,A);
            A = _mm256_fmadd_pd(L22,U2x,A);
            A = _mm256_fmadd_pd(L23,U3x,A); 
            _mm256_storeu_pd(&rowi2[j], A);

            A =  _mm256_loadu_pd( &rowi3[j] );
            A = _mm256_fmadd_pd(L30,U0x,A);
            A = _mm256_fmadd_pd(L31,U1x,A);
            A = _mm256_fmadd_pd(L32,U2x,A);
            A = _mm256_fmadd_pd(L33,U3x,A); 
            _mm256_storeu_pd(&rowi3[j], A);                                 
        }
    }


    {
        double U00 = rowk0[N-2];
        double U10 = rowk1[N-2];
        double U20 = rowk2[N-2];
        double U30 = rowk3[N-2];
        __m256d U00v = _mm256_setr_pd(U00,U10,U20,U30);

        double U01 = rowk0[N-1];
        double U11 = rowk1[N-1];
        double U21 = rowk2[N-1];
        double U31 = rowk3[N-1];
        __m256d U01v = _mm256_setr_pd(U01,U11,U21,U31);       
        for(unsigned int i = k1; i < this->N-2; i=i+4)
        {
            double* __restrict__ rowi0 = &v_[(i+0)*alignN];
            double* __restrict__ rowi1 = &v_[(i+1)*alignN];
            double* __restrict__ rowi2 = &v_[(i+2)*alignN];
            double* __restrict__ rowi3 = &v_[(i+3)*alignN];

            __m256d L0v = _mm256_loadu_pd(&rowi0[k0+0]);
            __m256d L1v = _mm256_loadu_pd(&rowi1[k0+0]);   
            __m256d L2v = _mm256_loadu_pd(&rowi2[k0+0]);
            __m256d L3v = _mm256_loadu_pd(&rowi3[k0+0]);

            __m256d UL0 = _mm256_mul_pd(U00v,L0v);
            __m256d UL1 = _mm256_mul_pd(U01v,L0v);
            __m256d UL2 = _mm256_mul_pd(U00v,L1v);
            __m256d UL3 = _mm256_mul_pd(U01v,L1v);
            transpose4x4_pd(UL0,UL1,UL2,UL3);
            __m256d result = _mm256_add_pd(UL0,UL1);
            result = _mm256_add_pd(result,UL2);
            result = _mm256_add_pd(result,UL3);
            __m128d r0 = _mm_loadu_pd(&rowi0[N-2]);
            __m128d r1 = _mm_loadu_pd(&rowi1[N-2]);
            r0 = _mm_sub_pd(r0,get_low(result));
            r1 = _mm_sub_pd(r1,get_high(result));
            _mm_storeu_pd(&rowi0[N-2],r0);
            _mm_storeu_pd(&rowi1[N-2],r1);


            __m256d UL4 = _mm256_mul_pd(U00v,L2v);
            __m256d UL5 = _mm256_mul_pd(U01v,L2v);
            __m256d UL6 = _mm256_mul_pd(U00v,L3v);
            __m256d UL7 = _mm256_mul_pd(U01v,L3v);
            transpose4x4_pd(UL4,UL5,UL6,UL7);
            __m256d result1 = _mm256_add_pd(UL4,UL5);
            result1 = _mm256_add_pd(result1,UL6);
            result1 = _mm256_add_pd(result1,UL7);
            __m128d r2 = _mm_loadu_pd(&rowi2[N-2]);
            __m128d r3 = _mm_loadu_pd(&rowi3[N-2]);
            r2 = _mm_sub_pd(r2,get_low(result1));
            r3 = _mm_sub_pd(r3,get_high(result1));
            _mm_storeu_pd(&rowi2[N-2],r2);
            _mm_storeu_pd(&rowi3[N-2],r3);
        }

        double* __restrict__ rowN1 = &v_[(N-1)*alignN];
        double* __restrict__ rowN2 = &v_[(N-2)*alignN];

        double L10 = rowN1[k0+0];
        double L11 = rowN1[k0+1];
        double L12 = rowN1[k0+2];
        double L13 = rowN1[k0+3];

        double L00 = rowN2[k0+0];
        double L01 = rowN2[k0+1];
        double L02 = rowN2[k0+2];
        double L03 = rowN2[k0+3];
        __m256d L10v = _mm256_set1_pd(-L10);
        __m256d L11v = _mm256_set1_pd(-L11);
        __m256d L12v = _mm256_set1_pd(-L12);
        __m256d L13v = _mm256_set1_pd(-L13);
        __m256d L00v = _mm256_set1_pd(-L00);
        __m256d L01v = _mm256_set1_pd(-L01);
        __m256d L02v = _mm256_set1_pd(-L02);
        __m256d L03v = _mm256_set1_pd(-L03);
        for(unsigned int i = k1; i < this->N-2; i=i+4)
        {
            __m256d rowk0v = _mm256_loadu_pd(&rowk0[i+0]);
            __m256d rowk1v = _mm256_loadu_pd(&rowk1[i+0]);
            __m256d rowk2v = _mm256_loadu_pd(&rowk2[i+0]);
            __m256d rowk3v = _mm256_loadu_pd(&rowk3[i+0]);

            __m256d rowN2v = _mm256_loadu_pd(&rowN2[i+0]);
            rowN2v = _mm256_fmadd_pd(rowk0v,L00v,rowN2v);
            rowN2v = _mm256_fmadd_pd(rowk1v,L01v,rowN2v);
            rowN2v = _mm256_fmadd_pd(rowk2v,L02v,rowN2v);
            rowN2v = _mm256_fmadd_pd(rowk3v,L03v,rowN2v);
            _mm256_storeu_pd(&rowN2[i+0],rowN2v);

            __m256d rowN1v = _mm256_loadu_pd(&rowN1[i+0]);
            rowN1v = _mm256_fmadd_pd(rowk0v,L10v,rowN1v);
            rowN1v = _mm256_fmadd_pd(rowk1v,L11v,rowN1v);
            rowN1v = _mm256_fmadd_pd(rowk2v,L12v,rowN1v);
            rowN1v = _mm256_fmadd_pd(rowk3v,L13v,rowN1v);
            _mm256_storeu_pd(&rowN1[i+0],rowN1v);           
        }
        {

            __m128d L10vv = _mm_set1_pd(-L10);
            __m128d L11vv = _mm_set1_pd(-L11);
            __m128d L12vv = _mm_set1_pd(-L12);
            __m128d L13vv = _mm_set1_pd(-L13);
            __m128d L00vv = _mm_set1_pd(-L00);
            __m128d L01vv = _mm_set1_pd(-L01);
            __m128d L02vv = _mm_set1_pd(-L02);
            __m128d L03vv = _mm_set1_pd(-L03);        
            unsigned int i = N-2;
            __m128d rowk0v = _mm_loadu_pd(&rowk0[i+0]);
            __m128d rowk1v = _mm_loadu_pd(&rowk1[i+0]);
            __m128d rowk2v = _mm_loadu_pd(&rowk2[i+0]);
            __m128d rowk3v = _mm_loadu_pd(&rowk3[i+0]);

            __m128d rowN2v = _mm_loadu_pd(&rowN2[i+0]);
            rowN2v = _mm_fmadd_pd(rowk0v,L00vv,rowN2v);
            rowN2v = _mm_fmadd_pd(rowk1v,L01vv,rowN2v);
            rowN2v = _mm_fmadd_pd(rowk2v,L02vv,rowN2v);
            rowN2v = _mm_fmadd_pd(rowk3v,L03vv,rowN2v);
            _mm_storeu_pd(&rowN2[i+0],rowN2v);

            __m128d rowN1v = _mm_loadu_pd(&rowN1[i+0]);
            rowN1v = _mm_fmadd_pd(rowk0v,L10vv,rowN1v);
            rowN1v = _mm_fmadd_pd(rowk1v,L11vv,rowN1v);
            rowN1v = _mm_fmadd_pd(rowk2v,L12vv,rowN1v);
            rowN1v = _mm_fmadd_pd(rowk3v,L13vv,rowN1v);
            _mm_storeu_pd(&rowN1[i+0],rowN1v);            
        }
    }

}

void FastChemistry::LUsolver::UpdateL22U22_Vec2_3
(
    unsigned int k0,
    unsigned int k1
)
{
    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];
    for(unsigned int i = k1; i < this->N-3; i=i+4)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];
        double* __restrict__ rowi3  = &v_[(i+3)*alignN];

        const double Li30 = rowi3[k0+0];
        const double Li31 = rowi3[k0+1];
        const double Li32 = rowi3[k0+2];
        const double Li33 = rowi3[k0+3];

        const double Li20 = rowi2[k0+0];
        const double Li21 = rowi2[k0+1];
        const double Li22 = rowi2[k0+2];
        const double Li23 = rowi2[k0+3];

        const double Li10 = rowi1[k0+0];
        const double Li11 = rowi1[k0+1];
        const double Li12 = rowi1[k0+2];
        const double Li13 = rowi1[k0+3];

        const double Li00 = rowi0[k0+0];
        const double Li01 = rowi0[k0+1];
        const double Li02 = rowi0[k0+2];
        const double Li03 = rowi0[k0+3];
        __m256d L00 = _mm256_set1_pd( -Li00);
        __m256d L01 = _mm256_set1_pd( -Li01);
        __m256d L02 = _mm256_set1_pd( -Li02);
        __m256d L03 = _mm256_set1_pd( -Li03);
        __m256d L10 = _mm256_set1_pd( -Li10);
        __m256d L11 = _mm256_set1_pd( -Li11);
        __m256d L12 = _mm256_set1_pd( -Li12);
        __m256d L13 = _mm256_set1_pd( -Li13);
        __m256d L20 = _mm256_set1_pd( -Li20);
        __m256d L21 = _mm256_set1_pd( -Li21);
        __m256d L22 = _mm256_set1_pd( -Li22);
        __m256d L23 = _mm256_set1_pd( -Li23);  
        __m256d L30 = _mm256_set1_pd( -Li30);
        __m256d L31 = _mm256_set1_pd( -Li31);
        __m256d L32 = _mm256_set1_pd( -Li32);
        __m256d L33 = _mm256_set1_pd( -Li33);
        for(unsigned int j = k1; j < this->N-3; j=j+4)
        {
            __m256d U0x = _mm256_loadu_pd(&rowk0[j]);
            __m256d U1x = _mm256_loadu_pd(&rowk1[j]);
            __m256d U2x = _mm256_loadu_pd(&rowk2[j]);
            __m256d U3x = _mm256_loadu_pd(&rowk3[j]);

            __m256d A = _mm256_loadu_pd( &rowi0[j] );
            A = _mm256_fmadd_pd(L00,U0x,A);
            A = _mm256_fmadd_pd(L01,U1x,A);
            A = _mm256_fmadd_pd(L02,U2x,A);
            A = _mm256_fmadd_pd(L03,U3x,A);  
            _mm256_storeu_pd(&rowi0[j], A);

            A = _mm256_loadu_pd( &rowi1[j] );
            A = _mm256_fmadd_pd(L10,U0x,A);
            A = _mm256_fmadd_pd(L11,U1x,A);
            A = _mm256_fmadd_pd(L12,U2x,A);
            A = _mm256_fmadd_pd(L13,U3x,A); 
            _mm256_storeu_pd(&rowi1[j], A);

            A = _mm256_loadu_pd( &rowi2[j] );
            A = _mm256_fmadd_pd(L20,U0x,A);
            A = _mm256_fmadd_pd(L21,U1x,A);
            A = _mm256_fmadd_pd(L22,U2x,A);
            A = _mm256_fmadd_pd(L23,U3x,A); 
            _mm256_storeu_pd(&rowi2[j], A);

            A = _mm256_loadu_pd( &rowi3[j] );
            A = _mm256_fmadd_pd(L30,U0x,A);
            A = _mm256_fmadd_pd(L31,U1x,A);
            A = _mm256_fmadd_pd(L32,U2x,A);
            A = _mm256_fmadd_pd(L33,U3x,A); 
            _mm256_storeu_pd(&rowi3[j], A);
        }
    }

    {
        double U00 = rowk0[N-3];
        double U10 = rowk1[N-3];
        double U20 = rowk2[N-3];
        double U30 = rowk3[N-3];
        __m256d U00v = _mm256_setr_pd(U00,U10,U20,U30);

        double U01 = rowk0[N-2];
        double U11 = rowk1[N-2];
        double U21 = rowk2[N-2];
        double U31 = rowk3[N-2];
        __m256d U01v = _mm256_setr_pd(U01,U11,U21,U31);

        double U02 = rowk0[N-1];
        double U12 = rowk1[N-1];
        double U22 = rowk2[N-1];
        double U32 = rowk3[N-1];
        __m256d U02v = _mm256_setr_pd(U02,U12,U22,U32);

        for(unsigned int i = k1; i < this->N-3; i=i+4)
        {
            double* __restrict__ rowi0  = &v_[(i+0)*alignN];
            double* __restrict__ rowi1  = &v_[(i+1)*alignN];
            double* __restrict__ rowi2  = &v_[(i+2)*alignN];
            double* __restrict__ rowi3  = &v_[(i+3)*alignN];

            __m256d L00v = _mm256_loadu_pd(&rowi0[k0+0]);
            __m256d L01v = _mm256_loadu_pd(&rowi1[k0+0]);
            __m256d L02v = _mm256_loadu_pd(&rowi2[k0+0]);
            __m256d L03v = _mm256_loadu_pd(&rowi3[k0+0]);

            __m256d L00vU00v = _mm256_mul_pd(L00v,U00v);
            __m256d L01vU00v = _mm256_mul_pd(L01v,U00v);
            __m256d L02vU00v = _mm256_mul_pd(L02v,U00v);
            __m256d L03vU00v = _mm256_mul_pd(L03v,U00v);

            __m256d L00vU01v = _mm256_mul_pd(L00v,U01v);
            __m256d L01vU01v = _mm256_mul_pd(L01v,U01v);
            __m256d L02vU01v = _mm256_mul_pd(L02v,U01v);
            __m256d L03vU01v = _mm256_mul_pd(L03v,U01v);

            __m256d L00vU02v = _mm256_mul_pd(L00v,U02v);
            __m256d L01vU02v = _mm256_mul_pd(L01v,U02v);
            __m256d L02vU02v = _mm256_mul_pd(L02v,U02v);
            __m256d L03vU02v = _mm256_mul_pd(L03v,U02v);

            rowi0[N-3] = rowi0[N-3] - hsum4(L00vU00v);
            rowi0[N-2] = rowi0[N-2] - hsum4(L00vU01v);
            rowi0[N-1] = rowi0[N-1] - hsum4(L00vU02v);

            rowi1[N-3] = rowi1[N-3] - hsum4(L01vU00v);
            rowi1[N-2] = rowi1[N-2] - hsum4(L01vU01v);
            rowi1[N-1] = rowi1[N-1] - hsum4(L01vU02v);

            rowi2[N-3] = rowi2[N-3] - hsum4(L02vU00v);
            rowi2[N-2] = rowi2[N-2] - hsum4(L02vU01v);
            rowi2[N-1] = rowi2[N-1] - hsum4(L02vU02v);

            rowi3[N-3] = rowi3[N-3] - hsum4(L03vU00v);
            rowi3[N-2] = rowi3[N-2] - hsum4(L03vU01v);
            rowi3[N-1] = rowi3[N-1] - hsum4(L03vU02v); 
        }

        double* __restrict__ rowN1  = &v_[(N-1)*alignN];
        double* __restrict__ rowN2  = &v_[(N-2)*alignN];
        double* __restrict__ rowN3  = &v_[(N-3)*alignN];
        double L20 = rowN1[k0+0];
        double L21 = rowN1[k0+1];
        double L22 = rowN1[k0+2];
        double L23 = rowN1[k0+3];
        __m256d L20v = _mm256_set1_pd(-L20);
        __m256d L21v = _mm256_set1_pd(-L21);
        __m256d L22v = _mm256_set1_pd(-L22);
        __m256d L23v = _mm256_set1_pd(-L23);
        double L10 = rowN2[k0+0];
        double L11 = rowN2[k0+1];
        double L12 = rowN2[k0+2];
        double L13 = rowN2[k0+3];
        __m256d L10v = _mm256_set1_pd(-L10);
        __m256d L11v = _mm256_set1_pd(-L11);
        __m256d L12v = _mm256_set1_pd(-L12);
        __m256d L13v = _mm256_set1_pd(-L13);
        double L00 = rowN3[k0+0];
        double L01 = rowN3[k0+1];
        double L02 = rowN3[k0+2];
        double L03 = rowN3[k0+3];
        __m256d L00v = _mm256_set1_pd(-L00);
        __m256d L01v = _mm256_set1_pd(-L01);
        __m256d L02v = _mm256_set1_pd(-L02);
        __m256d L03v = _mm256_set1_pd(-L03);
        for(unsigned int i = k1; i < this->N-3; i=i+4)
        {
            __m256d rowk0v = _mm256_loadu_pd(&rowk0[i+0]);
            __m256d rowk1v = _mm256_loadu_pd(&rowk1[i+0]);
            __m256d rowk2v = _mm256_loadu_pd(&rowk2[i+0]);
            __m256d rowk3v = _mm256_loadu_pd(&rowk3[i+0]);

            __m256d rowN3v = _mm256_loadu_pd(&rowN3[i+0]);
            rowN3v = _mm256_fmadd_pd(rowk0v,L00v,rowN3v);
            rowN3v = _mm256_fmadd_pd(rowk1v,L01v,rowN3v);
            rowN3v = _mm256_fmadd_pd(rowk2v,L02v,rowN3v);
            rowN3v = _mm256_fmadd_pd(rowk3v,L03v,rowN3v);
            _mm256_storeu_pd(&rowN3[i+0],rowN3v);

            __m256d rowN2v = _mm256_loadu_pd(&rowN2[i+0]);
            rowN2v = _mm256_fmadd_pd(rowk0v,L10v,rowN2v);
            rowN2v = _mm256_fmadd_pd(rowk1v,L11v,rowN2v);
            rowN2v = _mm256_fmadd_pd(rowk2v,L12v,rowN2v);
            rowN2v = _mm256_fmadd_pd(rowk3v,L13v,rowN2v);
            _mm256_storeu_pd(&rowN2[i+0],rowN2v);
  
            __m256d rowN1v = _mm256_loadu_pd(&rowN1[i+0]);
            rowN1v = _mm256_fmadd_pd(rowk0v,L20v,rowN1v);
            rowN1v = _mm256_fmadd_pd(rowk1v,L21v,rowN1v);
            rowN1v = _mm256_fmadd_pd(rowk2v,L22v,rowN1v);
            rowN1v = _mm256_fmadd_pd(rowk3v,L23v,rowN1v);
            _mm256_storeu_pd(&rowN1[i+0],rowN1v);

        }

        {
            __m128d L20vv = _mm_set1_pd(-L20);
            __m128d L21vv = _mm_set1_pd(-L21);
            __m128d L22vv = _mm_set1_pd(-L22);
            __m128d L23vv = _mm_set1_pd(-L23);
            __m128d L10vv = _mm_set1_pd(-L10);
            __m128d L11vv = _mm_set1_pd(-L11);
            __m128d L12vv = _mm_set1_pd(-L12);
            __m128d L13vv = _mm_set1_pd(-L13);
            __m128d L00vv = _mm_set1_pd(-L00);
            __m128d L01vv = _mm_set1_pd(-L01);
            __m128d L02vv = _mm_set1_pd(-L02);
            __m128d L03vv = _mm_set1_pd(-L03);        
            unsigned int i = N-3;
            __m128d rowk0v = _mm_loadu_pd(&rowk0[i+0]);
            __m128d rowk1v = _mm_loadu_pd(&rowk1[i+0]);
            __m128d rowk2v = _mm_loadu_pd(&rowk2[i+0]);
            __m128d rowk3v = _mm_loadu_pd(&rowk3[i+0]);

            __m128d rowN3v = _mm_loadu_pd(&rowN3[i+0]);
            rowN3v = _mm_fmadd_pd(rowk0v,L00vv,rowN3v);
            rowN3v = _mm_fmadd_pd(rowk1v,L01vv,rowN3v);
            rowN3v = _mm_fmadd_pd(rowk2v,L02vv,rowN3v);
            rowN3v = _mm_fmadd_pd(rowk3v,L03vv,rowN3v);
            _mm_storeu_pd(&rowN3[i+0],rowN3v);

            __m128d rowN2v = _mm_loadu_pd(&rowN2[i+0]);
            rowN2v = _mm_fmadd_pd(rowk0v,L10vv,rowN2v);
            rowN2v = _mm_fmadd_pd(rowk1v,L11vv,rowN2v);
            rowN2v = _mm_fmadd_pd(rowk2v,L12vv,rowN2v);
            rowN2v = _mm_fmadd_pd(rowk3v,L13vv,rowN2v);
            _mm_storeu_pd(&rowN2[i+0],rowN2v);

            __m128d rowN1v = _mm_loadu_pd(&rowN1[i+0]);
            rowN1v = _mm_fmadd_pd(rowk0v,L20vv,rowN1v);
            rowN1v = _mm_fmadd_pd(rowk1v,L21vv,rowN1v);
            rowN1v = _mm_fmadd_pd(rowk2v,L22vv,rowN1v);
            rowN1v = _mm_fmadd_pd(rowk3v,L23vv,rowN1v);
            _mm_storeu_pd(&rowN1[i+0],rowN1v);

            rowN3[i+2] = rowN3[i+2] - rowk0[i+2]*L00 - rowk1[i+2]*L01 - rowk2[i+2]*L02 - rowk3[i+2]*L03;
            rowN2[i+2] = rowN2[i+2] - rowk0[i+2]*L10 - rowk1[i+2]*L11 - rowk2[i+2]*L12 - rowk3[i+2]*L13;
            rowN1[i+2] = rowN1[i+2] - rowk0[i+2]*L20 - rowk1[i+2]*L21 - rowk2[i+2]*L22 - rowk3[i+2]*L23;
        }

    }
}

void FastChemistry::LUsolver::UpdateL22U22_Vec_38_00
(
    unsigned int k0,
    unsigned int k1,
    unsigned int rowTimes,
    unsigned int colTimes
)
{
    unsigned int rowEnd = rowTimes*3+k1;
    unsigned int colEnd = colTimes*8+k1;

    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];
    for(unsigned int i = k1; i<rowEnd; i=i+3)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];

        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = -load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);

        __m256d L203v = -load256d(&rowi2[k0+0]);
        __m256d L20 = _mm256_permute4x64_pd(L203v, 0x00);
        __m256d L21 = _mm256_permute4x64_pd(L203v, 0x55);
        __m256d L22 = _mm256_permute4x64_pd(L203v, 0xAA);
        __m256d L23 = _mm256_permute4x64_pd(L203v, 0xFF);
        for(unsigned int j=k1; j<colEnd; j=j+8)
        {
            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14
            __m256d C = load256d(&rowi2[j]);        //15

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fmadd_pd(L00,U003v,A);       //16
            B = _mm256_fmadd_pd(L10,U003v,B);       //16
            C = _mm256_fmadd_pd(L20,U003v,C);       //15

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fmadd_pd(L01,U103v,A);       //16
            B = _mm256_fmadd_pd(L11,U103v,B);       //16
            C = _mm256_fmadd_pd(L21,U103v,C);       //15

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fmadd_pd(L02,U203v,A);       //16
            B = _mm256_fmadd_pd(L12,U203v,B);       //16
            C = _mm256_fmadd_pd(L22,U203v,C);       //15

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14
            C = _mm256_fmadd_pd(L23,U303v,C);       //13
            store256d(&rowi2[j], C);                //12

            
            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d E = load256d(&rowi1[j+4]);      //14
            __m256d F = load256d(&rowi2[j+4]);      //15

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fmadd_pd(L00,U047v,D);       //11
            E = _mm256_fmadd_pd(L10,U047v,E);       //10
            F = _mm256_fmadd_pd(L20,U047v,F);       //10

            __m256d U147v = load256d(&rowk1[j+4]);  //11
            D = _mm256_fmadd_pd(L01,U147v,D);       //11
            E = _mm256_fmadd_pd(L11,U147v,E);       //10
            F = _mm256_fmadd_pd(L21,U147v,F);       //10

            __m256d U247v = load256d(&rowk2[j+4]);  //11
            D = _mm256_fmadd_pd(L02,U247v,D);       //10
            E = _mm256_fmadd_pd(L12,U247v,E);       //10
            F = _mm256_fmadd_pd(L22,U247v,F);       //10

            __m256d U347v = load256d(&rowk3[j+4]);  //11
            D = _mm256_fmadd_pd(L03,U347v,D);       //11
            store256d(&rowi0[j+4], D);              //10
            E = _mm256_fmadd_pd(L13,U347v,E);       //9
            store256d(&rowi1[j+4], E);              //8
            F = _mm256_fmadd_pd(L23,U347v,F);       //9
            store256d(&rowi2[j+4], F);              //8
        }
    }
}

void FastChemistry::LUsolver::UpdateL22U22_Vec_38_10
(
    unsigned int k0,
    unsigned int k1,
    unsigned int rowTimes,
    unsigned int colTimes
)
{
    unsigned int rowEnd = rowTimes*3+k1;
    unsigned int colEnd = colTimes*8+k1;

    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];
    for(unsigned int i = k1; i<rowEnd; i=i+3)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];

        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = -load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);

        __m256d L203v = -load256d(&rowi2[k0+0]);
        __m256d L20 = _mm256_permute4x64_pd(L203v, 0x00);
        __m256d L21 = _mm256_permute4x64_pd(L203v, 0x55);
        __m256d L22 = _mm256_permute4x64_pd(L203v, 0xAA);
        __m256d L23 = _mm256_permute4x64_pd(L203v, 0xFF);
        for(unsigned int j=k1; j<colEnd; j=j+8)
        {
            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14
            __m256d C = load256d(&rowi2[j]);        //15

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fmadd_pd(L00,U003v,A);       //16
            B = _mm256_fmadd_pd(L10,U003v,B);       //16
            C = _mm256_fmadd_pd(L20,U003v,C);       //15

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fmadd_pd(L01,U103v,A);       //16
            B = _mm256_fmadd_pd(L11,U103v,B);       //16
            C = _mm256_fmadd_pd(L21,U103v,C);       //15

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fmadd_pd(L02,U203v,A);       //16
            B = _mm256_fmadd_pd(L12,U203v,B);       //16
            C = _mm256_fmadd_pd(L22,U203v,C);       //15

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14
            C = _mm256_fmadd_pd(L23,U303v,C);       //13
            store256d(&rowi2[j], C);                //12

            
            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d E = load256d(&rowi1[j+4]);      //14
            __m256d F = load256d(&rowi2[j+4]);      //15

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fmadd_pd(L00,U047v,D);       //11
            E = _mm256_fmadd_pd(L10,U047v,E);       //10
            F = _mm256_fmadd_pd(L20,U047v,F);       //10

            __m256d U147v = load256d(&rowk1[j+4]);  //11
            D = _mm256_fmadd_pd(L01,U147v,D);       //11
            E = _mm256_fmadd_pd(L11,U147v,E);       //10
            F = _mm256_fmadd_pd(L21,U147v,F);       //10

            __m256d U247v = load256d(&rowk2[j+4]);  //11
            D = _mm256_fmadd_pd(L02,U247v,D);       //10
            E = _mm256_fmadd_pd(L12,U247v,E);       //10
            F = _mm256_fmadd_pd(L22,U247v,F);       //10

            __m256d U347v = load256d(&rowk3[j+4]);  //11
            D = _mm256_fmadd_pd(L03,U347v,D);       //11
            store256d(&rowi0[j+4], D);              //10
            E = _mm256_fmadd_pd(L13,U347v,E);       //9
            store256d(&rowi1[j+4], E);              //8
            F = _mm256_fmadd_pd(L23,U347v,F);       //9
            store256d(&rowi2[j+4], F);              //8
        }
    }
    {
        unsigned int i=rowEnd;
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];

        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        for(unsigned int j=k1; j<colEnd; j=j+8)
        {
            __m256d A = load256d(&rowi0[j]);        //13
            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fmadd_pd(L00,U003v,A);       //16
            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fmadd_pd(L01,U103v,A);       //16
            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fmadd_pd(L02,U203v,A);       //16
            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15

            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fmadd_pd(L00,U047v,D);       //11
            __m256d U147v = load256d(&rowk1[j+4]);  //11
            D = _mm256_fmadd_pd(L01,U147v,D);       //11
            __m256d U247v = load256d(&rowk2[j+4]);  //11
            D = _mm256_fmadd_pd(L02,U247v,D);       //10
            __m256d U347v = load256d(&rowk3[j+4]);  //11
            D = _mm256_fmadd_pd(L03,U347v,D);       //11
            store256d(&rowi0[j+4], D);              //10
        }
    }
}

void FastChemistry::LUsolver::UpdateL22U22_Vec_38_20
(
    unsigned int k0,
    unsigned int k1,
    unsigned int rowTimes,
    unsigned int colTimes
)
{
    unsigned int rowEnd = rowTimes*3+k1;
    unsigned int colEnd = colTimes*8+k1;

    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];
    for(unsigned int i = k1; i<rowEnd; i=i+3)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];

        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = -load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);

        __m256d L203v = -load256d(&rowi2[k0+0]);
        __m256d L20 = _mm256_permute4x64_pd(L203v, 0x00);
        __m256d L21 = _mm256_permute4x64_pd(L203v, 0x55);
        __m256d L22 = _mm256_permute4x64_pd(L203v, 0xAA);
        __m256d L23 = _mm256_permute4x64_pd(L203v, 0xFF);
        for(unsigned int j=k1; j<colEnd; j=j+8)
        {

            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14
            __m256d C = load256d(&rowi2[j]);        //15

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fmadd_pd(L00,U003v,A);       //16
            B = _mm256_fmadd_pd(L10,U003v,B);       //16
            C = _mm256_fmadd_pd(L20,U003v,C);       //15

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fmadd_pd(L01,U103v,A);       //16
            B = _mm256_fmadd_pd(L11,U103v,B);       //16
            C = _mm256_fmadd_pd(L21,U103v,C);       //15

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fmadd_pd(L02,U203v,A);       //16
            B = _mm256_fmadd_pd(L12,U203v,B);       //16
            C = _mm256_fmadd_pd(L22,U203v,C);       //15

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14
            C = _mm256_fmadd_pd(L23,U303v,C);       //13
            store256d(&rowi2[j], C);                //12

            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d E = load256d(&rowi1[j+4]);      //14
            __m256d F = load256d(&rowi2[j+4]);      //15

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fmadd_pd(L00,U047v,D);       //11
            E = _mm256_fmadd_pd(L10,U047v,E);       //10
            F = _mm256_fmadd_pd(L20,U047v,F);       //10

            __m256d U147v = load256d(&rowk1[j+4]);  //11
            D = _mm256_fmadd_pd(L01,U147v,D);       //11
            E = _mm256_fmadd_pd(L11,U147v,E);       //10
            F = _mm256_fmadd_pd(L21,U147v,F);       //10

            __m256d U247v = load256d(&rowk2[j+4]);  //11
            D = _mm256_fmadd_pd(L02,U247v,D);       //10
            E = _mm256_fmadd_pd(L12,U247v,E);       //10
            F = _mm256_fmadd_pd(L22,U247v,F);       //10

            __m256d U347v = load256d(&rowk3[j+4]);  //11
            D = _mm256_fmadd_pd(L03,U347v,D);       //11
            store256d(&rowi0[j+4], D);              //10
            E = _mm256_fmadd_pd(L13,U347v,E);       //9
            store256d(&rowi1[j+4], E);              //8
            F = _mm256_fmadd_pd(L23,U347v,F);       //9
            store256d(&rowi2[j+4], F);              //8
        }
    }
    {
        unsigned int i = rowEnd;
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];

        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = -load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);

        for(unsigned int j=k1; j<colEnd; j=j+8)
        {
            __m256d A = load256d(&rowi0[j]);
            __m256d B = load256d(&rowi1[j]);

            __m256d U003v = load256d(&rowk0[j]);
            A = _mm256_fmadd_pd(L00,U003v,A);
            B = _mm256_fmadd_pd(L10,U003v,B);

            __m256d U103v = load256d(&rowk1[j]);    
            A = _mm256_fmadd_pd(L01,U103v,A);
            B = _mm256_fmadd_pd(L11,U103v,B);

            __m256d U203v = load256d(&rowk2[j]);
            A = _mm256_fmadd_pd(L02,U203v,A);
            B = _mm256_fmadd_pd(L12,U203v,B);

            __m256d U303v = load256d(&rowk3[j]);    
            A = _mm256_fmadd_pd(L03,U303v,A);       
            store256d(&rowi0[j], A);                
            B = _mm256_fmadd_pd(L13,U303v,B);       
            store256d(&rowi1[j], B);                
 
            __m256d D = load256d(&rowi0[j+4]);      
            __m256d E = load256d(&rowi1[j+4]);      

            __m256d U047v = load256d(&rowk0[j+4]);  
            D = _mm256_fmadd_pd(L00,U047v,D);       
            E = _mm256_fmadd_pd(L10,U047v,E);       

            __m256d U147v = load256d(&rowk1[j+4]);  
            D = _mm256_fmadd_pd(L01,U147v,D);       
            E = _mm256_fmadd_pd(L11,U147v,E);       

            __m256d U247v = load256d(&rowk2[j+4]);  
            D = _mm256_fmadd_pd(L02,U247v,D);       
            E = _mm256_fmadd_pd(L12,U247v,E);       

            __m256d U347v = load256d(&rowk3[j+4]);  
            D = _mm256_fmadd_pd(L03,U347v,D);       
            store256d(&rowi0[j+4], D);              
            E = _mm256_fmadd_pd(L13,U347v,E);       
            store256d(&rowi1[j+4], E);              
        }
    }
}

void FastChemistry::LUsolver::UpdateL22U22_Vec_38_04
(
    unsigned int k0,
    unsigned int k1,
    unsigned int rowTimes,
    unsigned int colTimes
)
{
    unsigned int rowEnd = rowTimes*3+k1;
    unsigned int colEnd = colTimes*8+k1;

    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];
    for(unsigned int i = k1; i<rowEnd; i=i+3)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];

        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = -load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);

        __m256d L203v = -load256d(&rowi2[k0+0]);
        __m256d L20 = _mm256_permute4x64_pd(L203v, 0x00);
        __m256d L21 = _mm256_permute4x64_pd(L203v, 0x55);
        __m256d L22 = _mm256_permute4x64_pd(L203v, 0xAA);
        __m256d L23 = _mm256_permute4x64_pd(L203v, 0xFF);
        
        for(unsigned int j=k1; j<colEnd; j=j+8)
        {

            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14
            __m256d C = load256d(&rowi2[j]);        //15

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fmadd_pd(L00,U003v,A);       //16
            B = _mm256_fmadd_pd(L10,U003v,B);       //16
            C = _mm256_fmadd_pd(L20,U003v,C);       //15

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fmadd_pd(L01,U103v,A);       //16
            B = _mm256_fmadd_pd(L11,U103v,B);       //16
            C = _mm256_fmadd_pd(L21,U103v,C);       //15

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fmadd_pd(L02,U203v,A);       //16
            B = _mm256_fmadd_pd(L12,U203v,B);       //16
            C = _mm256_fmadd_pd(L22,U203v,C);       //15

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14
            C = _mm256_fmadd_pd(L23,U303v,C);       //13
            store256d(&rowi2[j], C);                //12
            
            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d E = load256d(&rowi1[j+4]);      //14
            __m256d F = load256d(&rowi2[j+4]);      //15

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fmadd_pd(L00,U047v,D);       //11
            E = _mm256_fmadd_pd(L10,U047v,E);       //10
            F = _mm256_fmadd_pd(L20,U047v,F);       //10

            __m256d U147v = load256d(&rowk1[j+4]);  //11
            D = _mm256_fmadd_pd(L01,U147v,D);       //11
            E = _mm256_fmadd_pd(L11,U147v,E);       //10
            F = _mm256_fmadd_pd(L21,U147v,F);       //10

            __m256d U247v = load256d(&rowk2[j+4]);  //11
            D = _mm256_fmadd_pd(L02,U247v,D);       //10
            E = _mm256_fmadd_pd(L12,U247v,E);       //10
            F = _mm256_fmadd_pd(L22,U247v,F);       //10

            __m256d U347v = load256d(&rowk3[j+4]);  //11
            D = _mm256_fmadd_pd(L03,U347v,D);       //11
            store256d(&rowi0[j+4], D);              //10
            E = _mm256_fmadd_pd(L13,U347v,E);       //9
            store256d(&rowi1[j+4], E);              //8
            F = _mm256_fmadd_pd(L23,U347v,F);       //9
            store256d(&rowi2[j+4], F);              //8
        }
        {
            unsigned int j = colEnd;
            __m256d A = load256d(&rowi0[j]);        
            __m256d B = load256d(&rowi1[j]);        
            __m256d C = load256d(&rowi2[j]);        

            __m256d U003v = load256d(&rowk0[j]);    
            A = _mm256_fmadd_pd(L00,U003v,A);       
            B = _mm256_fmadd_pd(L10,U003v,B);       
            C = _mm256_fmadd_pd(L20,U003v,C);       

            __m256d U103v = load256d(&rowk1[j]);    
            A = _mm256_fmadd_pd(L01,U103v,A);       
            B = _mm256_fmadd_pd(L11,U103v,B);       
            C = _mm256_fmadd_pd(L21,U103v,C);       

            __m256d U203v = load256d(&rowk2[j]);    
            A = _mm256_fmadd_pd(L02,U203v,A);       
            B = _mm256_fmadd_pd(L12,U203v,B);       
            C = _mm256_fmadd_pd(L22,U203v,C);       

            __m256d U303v = load256d(&rowk3[j]);
            A = _mm256_fmadd_pd(L03,U303v,A);
            B = _mm256_fmadd_pd(L13,U303v,B);
            C = _mm256_fmadd_pd(L23,U303v,C);
            store256d(&rowi0[j], A);
            store256d(&rowi1[j], B);
            store256d(&rowi2[j], C);       
        }
    }
}

void FastChemistry::LUsolver::UpdateL22U22_Vec_38_14
(
    unsigned int k0,
    unsigned int k1,
    unsigned int rowTimes,
    unsigned int colTimes
)
{
    unsigned int rowEnd = rowTimes*3+k1;
    unsigned int colEnd = colTimes*8+k1;

    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];
    for(unsigned int i = k1; i<rowEnd; i=i+3)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];

        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = -load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);

        __m256d L203v = -load256d(&rowi2[k0+0]);
        __m256d L20 = _mm256_permute4x64_pd(L203v, 0x00);
        __m256d L21 = _mm256_permute4x64_pd(L203v, 0x55);
        __m256d L22 = _mm256_permute4x64_pd(L203v, 0xAA);
        __m256d L23 = _mm256_permute4x64_pd(L203v, 0xFF);
        
        for(unsigned int j=k1; j<colEnd; j=j+8)
        {

            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14
            __m256d C = load256d(&rowi2[j]);        //15

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fmadd_pd(L00,U003v,A);       //16
            B = _mm256_fmadd_pd(L10,U003v,B);       //16
            C = _mm256_fmadd_pd(L20,U003v,C);       //15

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fmadd_pd(L01,U103v,A);       //16
            B = _mm256_fmadd_pd(L11,U103v,B);       //16
            C = _mm256_fmadd_pd(L21,U103v,C);       //15

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fmadd_pd(L02,U203v,A);       //16
            B = _mm256_fmadd_pd(L12,U203v,B);       //16
            C = _mm256_fmadd_pd(L22,U203v,C);       //15

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14
            C = _mm256_fmadd_pd(L23,U303v,C);       //13
            store256d(&rowi2[j], C);                //12

            
            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d E = load256d(&rowi1[j+4]);      //14
            __m256d F = load256d(&rowi2[j+4]);      //15

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fmadd_pd(L00,U047v,D);       //11
            E = _mm256_fmadd_pd(L10,U047v,E);       //10
            F = _mm256_fmadd_pd(L20,U047v,F);       //10

            __m256d U147v = load256d(&rowk1[j+4]);  //11
            D = _mm256_fmadd_pd(L01,U147v,D);       //11
            E = _mm256_fmadd_pd(L11,U147v,E);       //10
            F = _mm256_fmadd_pd(L21,U147v,F);       //10

            __m256d U247v = load256d(&rowk2[j+4]);  //11
            D = _mm256_fmadd_pd(L02,U247v,D);       //10
            E = _mm256_fmadd_pd(L12,U247v,E);       //10
            F = _mm256_fmadd_pd(L22,U247v,F);       //10

            __m256d U347v = load256d(&rowk3[j+4]);  //11
            D = _mm256_fmadd_pd(L03,U347v,D);       //11
            store256d(&rowi0[j+4], D);              //10
            E = _mm256_fmadd_pd(L13,U347v,E);       //9
            store256d(&rowi1[j+4], E);              //8
            F = _mm256_fmadd_pd(L23,U347v,F);       //9
            store256d(&rowi2[j+4], F);              //8
        }
        {
            unsigned int j=colEnd;
            __m256d A = load256d(&rowi0[j]);        
            __m256d B = load256d(&rowi1[j]);        
            __m256d C = load256d(&rowi2[j]);        

            __m256d U003v = load256d(&rowk0[j]);    
            A = _mm256_fmadd_pd(L00,U003v,A);       
            B = _mm256_fmadd_pd(L10,U003v,B);       
            C = _mm256_fmadd_pd(L20,U003v,C);       

            __m256d U103v = load256d(&rowk1[j]);    
            A = _mm256_fmadd_pd(L01,U103v,A);       
            B = _mm256_fmadd_pd(L11,U103v,B);       
            C = _mm256_fmadd_pd(L21,U103v,C);       

            __m256d U203v = load256d(&rowk2[j]);    
            A = _mm256_fmadd_pd(L02,U203v,A);       
            B = _mm256_fmadd_pd(L12,U203v,B);       
            C = _mm256_fmadd_pd(L22,U203v,C);       

            __m256d U303v = load256d(&rowk3[j]);
            A = _mm256_fmadd_pd(L03,U303v,A);
            B = _mm256_fmadd_pd(L13,U303v,B);
            C = _mm256_fmadd_pd(L23,U303v,C);
            store256d(&rowi0[j], A);
            store256d(&rowi1[j], B);
            store256d(&rowi2[j], C);       
        }
    }
    {
        unsigned int i = rowEnd;
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];

        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        for(unsigned int j=k1; j<colEnd; j=j+8)
        {

            __m256d A = load256d(&rowi0[j]);        
            __m256d U003v = load256d(&rowk0[j]);    
            A = _mm256_fmadd_pd(L00,U003v,A);       
            __m256d U103v = load256d(&rowk1[j]);    
            A = _mm256_fmadd_pd(L01,U103v,A);       
            __m256d U203v = load256d(&rowk2[j]);    
            A = _mm256_fmadd_pd(L02,U203v,A);       
            __m256d U303v = load256d(&rowk3[j]);    
            A = _mm256_fmadd_pd(L03,U303v,A);       
            store256d(&rowi0[j], A);                

            __m256d D = load256d(&rowi0[j+4]);      
            __m256d U047v = load256d(&rowk0[j+4]);  
            D = _mm256_fmadd_pd(L00,U047v,D);       
            __m256d U147v = load256d(&rowk1[j+4]);  
            D = _mm256_fmadd_pd(L01,U147v,D);       
            __m256d U247v = load256d(&rowk2[j+4]);  
            D = _mm256_fmadd_pd(L02,U247v,D);       
            __m256d U347v = load256d(&rowk3[j+4]);  
            D = _mm256_fmadd_pd(L03,U347v,D);       
            store256d(&rowi0[j+4], D);              
        }
        {
            unsigned int j=colEnd;
            __m256d A = load256d(&rowi0[j]);
            __m256d U003v = load256d(&rowk0[j]);    
            A = _mm256_fmadd_pd(L00,U003v,A);
            __m256d U103v = load256d(&rowk1[j]);    
            A = _mm256_fmadd_pd(L01,U103v,A);       
            __m256d U203v = load256d(&rowk2[j]);    
            A = _mm256_fmadd_pd(L02,U203v,A);       
            __m256d U303v = load256d(&rowk3[j]);
            A = _mm256_fmadd_pd(L03,U303v,A);
            store256d(&rowi0[j], A);
        }
    }
}

void FastChemistry::LUsolver::UpdateL22U22_Vec_38_24
(
    unsigned int k0,
    unsigned int k1,
    unsigned int rowTimes,
    unsigned int colTimes
)
{
    unsigned int rowEnd = rowTimes*3+k1;
    unsigned int colEnd = colTimes*8+k1;

    double* __restrict__ rowk0  = &v_[(k0+0)*alignN];
    double* __restrict__ rowk1  = &v_[(k0+1)*alignN];
    double* __restrict__ rowk2  = &v_[(k0+2)*alignN];
    double* __restrict__ rowk3  = &v_[(k0+3)*alignN];
    for(unsigned int i = k1; i<rowEnd; i=i+3)
    {
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];
        double* __restrict__ rowi2  = &v_[(i+2)*alignN];

        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = -load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);

        __m256d L203v = -load256d(&rowi2[k0+0]);
        __m256d L20 = _mm256_permute4x64_pd(L203v, 0x00);
        __m256d L21 = _mm256_permute4x64_pd(L203v, 0x55);
        __m256d L22 = _mm256_permute4x64_pd(L203v, 0xAA);
        __m256d L23 = _mm256_permute4x64_pd(L203v, 0xFF);
        
        for(unsigned int j=k1; j<colEnd; j=j+8)
        {

            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14
            __m256d C = load256d(&rowi2[j]);        //15

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fmadd_pd(L00,U003v,A);       //16
            B = _mm256_fmadd_pd(L10,U003v,B);       //16
            C = _mm256_fmadd_pd(L20,U003v,C);       //15

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fmadd_pd(L01,U103v,A);       //16
            B = _mm256_fmadd_pd(L11,U103v,B);       //16
            C = _mm256_fmadd_pd(L21,U103v,C);       //15

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fmadd_pd(L02,U203v,A);       //16
            B = _mm256_fmadd_pd(L12,U203v,B);       //16
            C = _mm256_fmadd_pd(L22,U203v,C);       //15

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14
            C = _mm256_fmadd_pd(L23,U303v,C);       //13
            store256d(&rowi2[j], C);                //12

            
            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d E = load256d(&rowi1[j+4]);      //14
            __m256d F = load256d(&rowi2[j+4]);      //15

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fmadd_pd(L00,U047v,D);       //11
            E = _mm256_fmadd_pd(L10,U047v,E);       //10
            F = _mm256_fmadd_pd(L20,U047v,F);       //10

            __m256d U147v = load256d(&rowk1[j+4]);  //11
            D = _mm256_fmadd_pd(L01,U147v,D);       //11
            E = _mm256_fmadd_pd(L11,U147v,E);       //10
            F = _mm256_fmadd_pd(L21,U147v,F);       //10

            __m256d U247v = load256d(&rowk2[j+4]);  //11
            D = _mm256_fmadd_pd(L02,U247v,D);       //10
            E = _mm256_fmadd_pd(L12,U247v,E);       //10
            F = _mm256_fmadd_pd(L22,U247v,F);       //10

            __m256d U347v = load256d(&rowk3[j+4]);  //11
            D = _mm256_fmadd_pd(L03,U347v,D);       //11
            store256d(&rowi0[j+4], D);              //10
            E = _mm256_fmadd_pd(L13,U347v,E);       //9
            store256d(&rowi1[j+4], E);              //8
            F = _mm256_fmadd_pd(L23,U347v,F);       //9
            store256d(&rowi2[j+4], F);              //8
        }
        {
            unsigned int j=colEnd;
            __m256d A = load256d(&rowi0[j]);        
            __m256d B = load256d(&rowi1[j]);        
            __m256d C = load256d(&rowi2[j]);        

            __m256d U003v = load256d(&rowk0[j]);    
            A = _mm256_fmadd_pd(L00,U003v,A);       
            B = _mm256_fmadd_pd(L10,U003v,B);       
            C = _mm256_fmadd_pd(L20,U003v,C);       

            __m256d U103v = load256d(&rowk1[j]);    
            A = _mm256_fmadd_pd(L01,U103v,A);       
            B = _mm256_fmadd_pd(L11,U103v,B);       
            C = _mm256_fmadd_pd(L21,U103v,C);       

            __m256d U203v = load256d(&rowk2[j]);    
            A = _mm256_fmadd_pd(L02,U203v,A);       
            B = _mm256_fmadd_pd(L12,U203v,B);       
            C = _mm256_fmadd_pd(L22,U203v,C);       

            __m256d U303v = load256d(&rowk3[j]);
            A = _mm256_fmadd_pd(L03,U303v,A);
            B = _mm256_fmadd_pd(L13,U303v,B);
            C = _mm256_fmadd_pd(L23,U303v,C);
            store256d(&rowi0[j], A);
            store256d(&rowi1[j], B);
            store256d(&rowi2[j], C);       
        }
    }
    {
        unsigned int i = rowEnd;
        double* __restrict__ rowi0  = &v_[(i+0)*alignN];
        double* __restrict__ rowi1  = &v_[(i+1)*alignN];


        __m256d L003v = -load256d(&rowi0[k0+0]);
        __m256d L00 = _mm256_permute4x64_pd(L003v, 0x00);
        __m256d L01 = _mm256_permute4x64_pd(L003v, 0x55);
        __m256d L02 = _mm256_permute4x64_pd(L003v, 0xAA);
        __m256d L03 = _mm256_permute4x64_pd(L003v, 0xFF);

        __m256d L103v = -load256d(&rowi1[k0+0]);
        __m256d L10 = _mm256_permute4x64_pd(L103v, 0x00);
        __m256d L11 = _mm256_permute4x64_pd(L103v, 0x55);
        __m256d L12 = _mm256_permute4x64_pd(L103v, 0xAA);
        __m256d L13 = _mm256_permute4x64_pd(L103v, 0xFF);
        
        for(unsigned int j=k1; j<colEnd; j=j+8)
        {

            __m256d A = load256d(&rowi0[j]);        //13
            __m256d B = load256d(&rowi1[j]);        //14

            __m256d U003v = load256d(&rowk0[j]);    //16
            A = _mm256_fmadd_pd(L00,U003v,A);       //16
            B = _mm256_fmadd_pd(L10,U003v,B);       //16

            __m256d U103v = load256d(&rowk1[j]);    //16
            A = _mm256_fmadd_pd(L01,U103v,A);       //16
            B = _mm256_fmadd_pd(L11,U103v,B);       //16

            __m256d U203v = load256d(&rowk2[j]);    //16
            A = _mm256_fmadd_pd(L02,U203v,A);       //16
            B = _mm256_fmadd_pd(L12,U203v,B);       //16

            __m256d U303v = load256d(&rowk3[j]);    //16
            A = _mm256_fmadd_pd(L03,U303v,A);       //16
            store256d(&rowi0[j], A);                //15
            B = _mm256_fmadd_pd(L13,U303v,B);       //15
            store256d(&rowi1[j], B);                //14

            
            __m256d D = load256d(&rowi0[j+4]);      //13
            __m256d E = load256d(&rowi1[j+4]);      //14

            __m256d U047v = load256d(&rowk0[j+4]);  //11
            D = _mm256_fmadd_pd(L00,U047v,D);       //11
            E = _mm256_fmadd_pd(L10,U047v,E);       //10

            __m256d U147v = load256d(&rowk1[j+4]);  //11
            D = _mm256_fmadd_pd(L01,U147v,D);       //11
            E = _mm256_fmadd_pd(L11,U147v,E);       //10

            __m256d U247v = load256d(&rowk2[j+4]);  //11
            D = _mm256_fmadd_pd(L02,U247v,D);       //10
            E = _mm256_fmadd_pd(L12,U247v,E);       //10

            __m256d U347v = load256d(&rowk3[j+4]);  //11
            D = _mm256_fmadd_pd(L03,U347v,D);       //11
            store256d(&rowi0[j+4], D);              //10
            E = _mm256_fmadd_pd(L13,U347v,E);       //9
            store256d(&rowi1[j+4], E);              //8
        }
        {
            unsigned int j=colEnd;
            __m256d A = load256d(&rowi0[j]);
            __m256d B = load256d(&rowi1[j]);

            __m256d U003v = load256d(&rowk0[j]);    
            A = _mm256_fmadd_pd(L00,U003v,A);       
            B = _mm256_fmadd_pd(L10,U003v,B);       

            __m256d U103v = load256d(&rowk1[j]);    
            A = _mm256_fmadd_pd(L01,U103v,A);       
            B = _mm256_fmadd_pd(L11,U103v,B);       

            __m256d U203v = load256d(&rowk2[j]);    
            A = _mm256_fmadd_pd(L02,U203v,A);       
            B = _mm256_fmadd_pd(L12,U203v,B);       

            __m256d U303v = load256d(&rowk3[j]);
            A = _mm256_fmadd_pd(L03,U303v,A);
            B = _mm256_fmadd_pd(L13,U303v,B);
            store256d(&rowi0[j], A);
            store256d(&rowi1[j], B);
        }
    }
}

void FastChemistry::LUsolver::permutation0
(
    unsigned int k0,
    unsigned int k1
)
{
    for(unsigned int j = 0; j<4;j++)
    {
        if(this->pivotIndice_[k0+j]==k0+j)
        {
            continue;
        }
        for(unsigned int i = 0; i <= this->N-4; i=i+4)
        {
            unsigned int iTarget = this->pivotIndice_[j+k0];
            double* __restrict__ rowiTarget = &v_[iTarget*alignN];
            double* __restrict__ rowj = &v_[(k0+j)*alignN];

            __m256d A = load256d(&rowiTarget[i]);
            __m256d B = load256d(&rowj[i]);
            store256d(&rowiTarget[i],B);
            store256d(&rowj[i],A);
        }
    }
}

void FastChemistry::LUsolver::permutation1
(
    unsigned int k0,
    unsigned int k1
)
{
    for(unsigned int j = 0; j<4;j++)
    {
        if(this->pivotIndice_[k0+j]==k0+j)
        {
            continue;
        }
        for(unsigned int i = 0; i <= this->N-4; i=i+4)
        {
            unsigned int iTarget = this->pivotIndice_[j+k0];
            double* __restrict__ rowiTarget = &v_[iTarget*alignN];
            double* __restrict__ rowj = &v_[(k0+j)*alignN];

            __m256d A = load256d(&rowiTarget[i]);
            __m256d B = load256d(&rowj[i]);
            store256d(&rowiTarget[i],B);
            store256d(&rowj[i],A);
        }
        {
            unsigned int i = this->N-1;
            unsigned int iTarget = this->pivotIndice_[j+k0];            
            double* __restrict__ rowiTarget = &v_[iTarget*alignN];
            double* __restrict__ rowj = &v_[(k0+j)*alignN];            
            double A = rowiTarget[i]; 
            rowiTarget[i] = rowj[i];
            rowj[i] = A;
        }
    }
}

void FastChemistry::LUsolver::permutation2
(
    unsigned int k0,
    unsigned int k1
)
{
    for(unsigned int j = 0; j<4;j++)
    {
        if(this->pivotIndice_[k0+j]==k0+j)
        {
            continue;
        }
        for(unsigned int i = 0; i <= this->N-4; i=i+4)
        {
            unsigned int iTarget = this->pivotIndice_[j+k0];
            double* __restrict__ rowiTarget = &v_[iTarget*alignN];
            double* __restrict__ rowj = &v_[(k0+j)*alignN];

            __m256d A = load256d(&rowiTarget[i]);
            __m256d B = load256d(&rowj[i]);
            store256d(&rowiTarget[i],B);
            store256d(&rowj[i],A);
        }
        {
            unsigned int i = this->N-2;
            unsigned int iTarget = this->pivotIndice_[j+k0];            
            double* __restrict__ rowiTarget = &v_[iTarget*alignN];
            double* __restrict__ rowj = &v_[(k0+j)*alignN];          

            double A0 = rowiTarget[i+0]; 
            double A1 = rowiTarget[i+1];
            rowiTarget[i+0] = rowj[i+0];
            rowiTarget[i+1] = rowj[i+1];
            rowj[i+0] = A0;
            rowj[i+1] = A1;
        }
    }
}

void FastChemistry::LUsolver::permutation3
(
    unsigned int k0,
    unsigned int k1
)
{
    for(unsigned int j = 0; j<4;j++)
    {
        if(this->pivotIndice_[k0+j]==k0+j)
        {
            continue;
        }
        for(unsigned int i = 0; i <= this->N-4; i=i+4)
        {
            unsigned int iTarget = this->pivotIndice_[j+k0];
            double* __restrict__ rowiTarget = &v_[iTarget*alignN];
            double* __restrict__ rowj = &v_[(k0+j)*alignN];

            __m256d A = load256d(&rowiTarget[i]);
            __m256d B = load256d(&rowj[i]);
            store256d(&rowiTarget[i],B);
            store256d(&rowj[i],A);
        }
        {
            unsigned int i = this->N-3;
            unsigned int iTarget = this->pivotIndice_[j+k0];            
            double* __restrict__ rowiTarget = &v_[iTarget*alignN];
            double* __restrict__ rowj = &v_[(k0+j)*alignN];          

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

void FastChemistry::LUsolver::permutation_2
(
)
{
    double Array0[4];    
    {
        if(this->pivotIndice_[N-2]==N-1)
        {
            for(unsigned int i = 0; i < this->N-2; i=i+4)
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

void FastChemistry::LUsolver::permutation_3
(
)
{
    double Array0[4];    
    {
        unsigned int k0 = this->N-3;
        for(unsigned int j = 0; j<3;j++)
        {
            if(this->pivotIndice_[k0+j]==k0+j)
            {
                continue;
            }

            for(unsigned int i = 0; i < k0; i=i+4)
            {
                unsigned int iTarget = this->pivotIndice_[j+k0];

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





void FastChemistry::LUsolver::xSolve
(
    double* __restrict__ b
)
{
    if(N<8)
    {
        this->xSolve_Serial(b);
        return;
    }
    if(Remain==0)
    {
        this->xSolve_Vec_0(b);
    }
    else if(Remain==1)
    {
        this->xSolve_Vec_1(b);
    }
    else if(Remain==2)
    {
        this->xSolve_Vec_2(b);
    }
    else
    {
        this->xSolve_Vec_3(b);
    }
}

void FastChemistry::LUsolver::xSolve_Vec_3
(
    double* __restrict__ b
)
{
    for(unsigned int j = 0; j<this->N;j++)
    {
        if(this->pivotIndice_[j]==j)
        {
            continue;
        }
        unsigned int jTarget = this->pivotIndice_[j];
        double temp = b[jTarget];
        b[jTarget] = b[j];
        b[j] = temp;
    }
    {
        double b0 = b[0];
        double b1 = b[1];
        double b2 = b[2];
        double b3 = b[3];

        double L10 = v_[1*alignN+0];
        double L20 = v_[2*alignN+0];
        double L21 = v_[2*alignN+1];
        double L30 = v_[3*alignN+0];
        double L31 = v_[3*alignN+1];
        double L32 = v_[3*alignN+2];
        b1 = std::fma(-L10, b0, b1);
        b2 = std::fma(-L20, b0, b2);
        b2 = std::fma(-L21, b1, b2);
        b3 = std::fma(-L30, b0, b3);
        b3 = std::fma(-L31, b1, b3);
        b3 = std::fma(-L32, b2, b3);
        b[1] = b1;
        b[2] = b2;
        b[3] = b3;
    }

    for(unsigned int i = 4; i < this->N-3; i=i+4)
    {
        unsigned int i0 = i+0;
        unsigned int i1 = i+1;
        unsigned int i2 = i+2;
        unsigned int i3 = i+3;

        __m256d sum = _mm256_loadu_pd(&b[i0]);
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ];
        double* __restrict__ Lrow2 = &v_[ i2 * alignN ];
        double* __restrict__ Lrow3 = &v_[ i3 * alignN ];

        for(unsigned int j = 0; j < i; j=j+4)
        {
            __m256d L0 = _mm256_loadu_pd(&Lrow0[j]);
            __m256d L1 = _mm256_loadu_pd(&Lrow1[j]);
            __m256d L2 = _mm256_loadu_pd(&Lrow2[j]);
            __m256d L3 = _mm256_loadu_pd(&Lrow3[j]);
            __m256d Y = _mm256_loadu_pd(&b[j+0]);
            __m256d L0Y = _mm256_mul_pd(L0,Y);
            __m256d L1Y = _mm256_mul_pd(L1,Y);
            __m256d L2Y = _mm256_mul_pd(L2,Y);
            __m256d L3Y = _mm256_mul_pd(L3,Y);
            __m256d temp1 = _mm256_hadd_pd(L0Y,L1Y);
            __m256d temp2 = _mm256_hadd_pd(L2Y,L3Y);
            __m256d temp1_shuffle = _mm256_permute4x64_pd(temp1,_MM_SHUFFLE(3,1,2,0));
            __m256d temp2_shuffle = _mm256_permute4x64_pd(temp2,_MM_SHUFFLE(3,1,2,0));
            __m256d temp3 = _mm256_hadd_pd(temp1_shuffle,temp2_shuffle);
            __m256d temp3_shuffle = _mm256_permute4x64_pd(temp3,_MM_SHUFFLE(3,1,2,0));
            sum = _mm256_add_pd(sum,-temp3_shuffle);
        }
        _mm256_storeu_pd(&b[i0],sum);
       
        double b0 = b[i0];
        double b1 = b[i1];
        double b2 = b[i2];
        double b3 = b[i3];
        double L10 = Lrow1[i0];
        double L20 = Lrow2[i0];
        double L21 = Lrow2[i1];
        double L30 = Lrow3[i0];
        double L31 = Lrow3[i1];
        double L32 = Lrow3[i2];
        b1 = std::fma(-L10, b0, b1);
        b2 = std::fma(-L20, b0, b2);
        b2 = std::fma(-L21, b1, b2);
        b3 = std::fma(-L30, b0, b3);
        b3 = std::fma(-L31, b1, b3);
        b3 = std::fma(-L32, b2, b3);
        b[i1] = b1;
        b[i2] = b2;
        b[i3] = b3;      
    }
    {
        unsigned int i0 = this->N-3;        
        unsigned int i1 = this->N-2;
        unsigned int i2 = this->N-1;

        __m256d sum = _mm256_set1_pd(0);

        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ];
        double* __restrict__ Lrow2 = &v_[ i2 * alignN ];
        for(unsigned int j = 0; j < i0; j=j+4)
        {
            __m256d L0 = _mm256_loadu_pd(&Lrow0[j]);
            __m256d L1 = _mm256_loadu_pd(&Lrow1[j]);
            __m256d L2 = _mm256_loadu_pd(&Lrow2[j]);
            __m256d Y = _mm256_loadu_pd(&b[j+0]);
            __m256d L0Y = _mm256_mul_pd(L0,Y);
            __m256d L1Y = _mm256_mul_pd(L1,Y);
            __m256d L2Y = _mm256_mul_pd(L2,Y);
            __m256d temp1 = _mm256_hadd_pd(L0Y,L1Y);
            __m256d temp2 = _mm256_hadd_pd(L2Y,L2Y);
            __m256d temp1_shuffle = _mm256_permute4x64_pd(temp1,_MM_SHUFFLE(3,1,2,0));
            __m256d temp2_shuffle = _mm256_permute4x64_pd(temp2,_MM_SHUFFLE(3,1,2,0));
            __m256d temp3 = _mm256_hadd_pd(temp1_shuffle,temp2_shuffle);
            __m256d temp3_shuffle = _mm256_permute4x64_pd(temp3,_MM_SHUFFLE(3,1,2,0));
            sum = _mm256_add_pd(sum,temp3_shuffle);
        }

        __m128d lo = _mm256_castpd256_pd128(sum);      
        __m128d hi = _mm256_extractf128_pd(sum, 1);    
        double s0 = _mm_cvtsd_f64(lo);  
        __m128d lo_hi = _mm_unpackhi_pd(lo, lo);  
        double s1 = _mm_cvtsd_f64(lo_hi);
        double s2 = _mm_cvtsd_f64(hi);

        b[i0] = b[i0] - s0;
        b[i1] = b[i1] - s1 - Lrow1[i0]*b[i0];
        b[i2] = b[i2] - s2 - Lrow2[i0]*b[i0] - Lrow2[i1]*b[i1];        
    }
    {

        double b3 = b[N-1];
        double b2 = b[N-2];
        double b1 = b[N-3];
        double b0 = b[N-4];

        double* __restrict__ rowN2 = &v_[(N-2)*alignN];
        double* __restrict__ rowN3 = &v_[(N-3)*alignN];
        double* __restrict__ rowN4 = &v_[(N-4)*alignN];
        const double invD33 = this->invD[N-1];
        const double invD22 = this->invD[N-2];
        const double invD11 = this->invD[N-3];
        const double invD00 = this->invD[N-4];
        double L23 = rowN2[N-1];
        double L12 = rowN3[N-2];        
        double L13 = rowN3[N-1];
        double L03 = rowN4[N-1];
        double L02 = rowN4[N-2];
        double L01 = rowN4[N-3];
        b3 = b3 * invD33;
        b2 = std::fma(-L23, b3, b2);
        b2 = b2 * invD22;
        b1 = std::fma(-L13, b3, b1);
        b1 = std::fma(-L12, b2, b1);
        b1 = b1 * invD11;
        b0 = std::fma(-L03, b3, b0);
        b0 = std::fma(-L02, b2, b0);
        b0 = std::fma(-L01, b1, b0);
        b0 = b0 * invD00;
        b[N-1] = b3;
        b[N-2] = b2;
        b[N-3] = b1;
        b[N-4] = b0;
    }
    for(int i = N-1-4; i >= 3; i=i-4)
    {
        const unsigned int i0 = i-3;
        const unsigned int i1 = i-2;
        const unsigned int i2 = i-1;
        const unsigned int i3 = i-0;

        double sum0 = b[i0];
        double sum1 = b[i1];
        double sum2 = b[i2];
        double sum3 = b[i3];
        __m256d result0 = _mm256_setzero_pd();
        __m256d result1 = _mm256_setzero_pd();
        __m256d result2 = _mm256_setzero_pd();
        __m256d result3 = _mm256_setzero_pd();
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ];
        double* __restrict__ Lrow2 = &v_[ i2 * alignN ];
        double* __restrict__ Lrow3 = &v_[ i3 * alignN ];

        for(unsigned int j = i+1; j<this->N; j=j+4)
        {
            unsigned int j0 = j + 0;
            __m256d bv = _mm256_loadu_pd(&b[j0]);
            __m256d A00v = _mm256_loadu_pd(&Lrow0[j0]);
            __m256d A10v = _mm256_loadu_pd(&Lrow1[j0]);
            __m256d A20v = _mm256_loadu_pd(&Lrow2[j0]);
            __m256d A30v = _mm256_loadu_pd(&Lrow3[j0]);

            result0 = _mm256_fmadd_pd(A00v,bv,result0);
            result1 = _mm256_fmadd_pd(A10v,bv,result1);
            result2 = _mm256_fmadd_pd(A20v,bv,result2);
            result3 = _mm256_fmadd_pd(A30v,bv,result3);        
        }              
        sum0 = sum0 - this->hsum4(result0);
        sum1 = sum1 - this->hsum4(result1);
        sum2 = sum2 - this->hsum4(result2);
        sum3 = sum3 - this->hsum4(result3); 
        const double invD33 = this->invD[i3];
        const double invD22 = this->invD[i2];
        const double invD11 = this->invD[i1];
        const double invD00 = this->invD[i0];
        double L23 = Lrow2[i3];
        double L12 = Lrow1[i2];
        double L13 = Lrow1[i3];
        double L01 = Lrow0[i1];
        double L02 = Lrow0[i2];
        double L03 = Lrow0[i3];

        double b3 = b[i3];
        double b2 = b[i2];
        double b1 = b[i1];
        double b0 = b[i0];

        b3 = sum3 * invD33;
        b2 = std::fma(-L23,b3,sum2);
        b2 = b2 * invD22;
        b1 = std::fma(-L12,b2,sum1);
        b1 = std::fma(-L13,b3,b1);
        b1 = b1 * invD11;
        b0 = std::fma(-L01,b1,sum0);
        b0 = std::fma(-L02,b2,b0);
        b0 = std::fma(-L03,b3,b0);
        b0 = b0 * invD00;

        b[i3] = b3;
        b[i2] = b2;
        b[i1] = b1;
        b[i0] = b0;    
    }
    {
        const unsigned int i0 = 0;
        const unsigned int i1 = 1;       
        const unsigned int i2 = 2;  
        double sum0 = b[i0];
        double sum1 = b[i1];
        double sum2 = b[i2];
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ];
        double* __restrict__ Lrow2 = &v_[ i2 * alignN ];
        __m256d result0 = _mm256_setzero_pd();
        __m256d result1 = _mm256_setzero_pd();
        __m256d result2 = _mm256_setzero_pd();

        for(unsigned int j = i2+1; j<this->N; j=j+4)
        {
            unsigned int j0 = j + 0;
            
            __m256d bv = _mm256_loadu_pd(&b[j0]);
            __m256d A0 = _mm256_loadu_pd(&Lrow0[j0]);
            __m256d A1 = _mm256_loadu_pd(&Lrow1[j0]);
            __m256d A2 = _mm256_loadu_pd(&Lrow2[j0]);
            result0 = _mm256_fmadd_pd(A0,bv,result0);
            result1 = _mm256_fmadd_pd(A1,bv,result1);
            result2 = _mm256_fmadd_pd(A2,bv,result2);
        }
        double r0 = hsum4(result0);
        double r1 = hsum4(result1);
        double r2 = hsum4(result2);
        sum0 = sum0 - r0;
        sum1 = sum1 - r1;
        sum2 = sum2 - r2;
        const double invD22 = this->invD[i2];
        const double invD11 = this->invD[i1];
        const double invD00 = this->invD[i0];
        b[i2] = (sum2) * invD22;          
        b[i1] = (sum1 - Lrow1[i2]*b[i2]) * invD11;   
        b[i0] = (sum0 - Lrow0[i2]*b[i2] - Lrow0[i1]*b[i1]) * invD00;   
    }
}
void FastChemistry::LUsolver::xSolve_Vec_2
(
    double* __restrict__ b
)
{
    for(unsigned int j = 0; j<this->N;j++)
    {
        if(this->pivotIndice_[j]==j)
        {
            continue;
        }

        unsigned int jTarget = this->pivotIndice_[j];

        double temp = b[jTarget];
        b[jTarget] = b[j];
        b[j] = temp;
    }
    {
        double b0 = b[0];
        double b1 = b[1];
        double b2 = b[2];
        double b3 = b[3];
        double L10 = v_[1*alignN+0];
        double L20 = v_[2*alignN+0];
        double L21 = v_[2*alignN+1];
        double L30 = v_[3*alignN+0];
        double L31 = v_[3*alignN+1];
        double L32 = v_[3*alignN+2];
        b1 = std::fma(-L10, b0, b1);
        b2 = std::fma(-L20, b0, b2);
        b2 = std::fma(-L21, b1, b2);
        b3 = std::fma(-L30, b0, b3);
        b3 = std::fma(-L31, b1, b3);
        b3 = std::fma(-L32, b2, b3);
        b[1] = b1;
        b[2] = b2;
        b[3] = b3;
    }
    for(unsigned int i = 4; i < this->N-2; i=i+4)
    {
        unsigned int i0 = i+0;
        unsigned int i1 = i+1;
        unsigned int i2 = i+2;
        unsigned int i3 = i+3;

        __m256d sum = _mm256_loadu_pd(&b[i0]);
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ];
        double* __restrict__ Lrow2 = &v_[ i2 * alignN ];
        double* __restrict__ Lrow3 = &v_[ i3 * alignN ];        
        for(unsigned int j = 0; j < i; j=j+4)
        {
            __m256d L0 = _mm256_loadu_pd(&Lrow0[j]);
            __m256d L1 = _mm256_loadu_pd(&Lrow1[j]);
            __m256d L2 = _mm256_loadu_pd(&Lrow2[j]);
            __m256d L3 = _mm256_loadu_pd(&Lrow3[j]);
            __m256d Y = _mm256_loadu_pd(&b[j+0]);
            __m256d L0Y = _mm256_mul_pd(L0,Y);
            __m256d L1Y = _mm256_mul_pd(L1,Y);
            __m256d L2Y = _mm256_mul_pd(L2,Y);
            __m256d L3Y = _mm256_mul_pd(L3,Y);
            
            __m256d temp1 = _mm256_hadd_pd(L0Y,L1Y);
            __m256d temp2 = _mm256_hadd_pd(L2Y,L3Y);
            __m256d temp1_shuffle = _mm256_permute4x64_pd(temp1,_MM_SHUFFLE(3,1,2,0));
            __m256d temp2_shuffle = _mm256_permute4x64_pd(temp2,_MM_SHUFFLE(3,1,2,0));
            __m256d temp3 = _mm256_hadd_pd(temp1_shuffle,temp2_shuffle);
            __m256d temp3_shuffle = _mm256_permute4x64_pd(temp3,_MM_SHUFFLE(3,1,2,0));
            sum = _mm256_add_pd(sum,-temp3_shuffle);
        }
        _mm256_storeu_pd(&b[i0],sum);    
        double b0 = b[i0];
        double b1 = b[i1];
        double b2 = b[i2];
        double b3 = b[i3];
        double L10 = Lrow1[i0];
        double L20 = Lrow2[i0];
        double L21 = Lrow2[i1];
        double L30 = Lrow3[i0];
        double L31 = Lrow3[i1];
        double L32 = Lrow3[i2];
        b1 = std::fma(-L10, b0, b1);
        b2 = std::fma(-L20, b0, b2);
        b2 = std::fma(-L21, b1, b2);
        b3 = std::fma(-L30, b0, b3);
        b3 = std::fma(-L31, b1, b3);
        b3 = std::fma(-L32, b2, b3);
        b[i1] = b1;
        b[i2] = b2;
        b[i3] = b3;    
    }

    {
        unsigned int i0 = this->N-2;
        unsigned int i1 = this->N-1;
        int remain = i0%4;
        __m256d sum1 = _mm256_set1_pd(0);
        __m256d sum2 = _mm256_set1_pd(0);
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ];       
        for(unsigned int j = 0; j < i0-remain; j=j+4)
        {
            __m256d L0 = _mm256_loadu_pd(&Lrow0[j]);
            __m256d L1 = _mm256_loadu_pd(&Lrow1[j]);
            __m256d Y = _mm256_loadu_pd(&b[j+0]);
            sum1 = _mm256_fmadd_pd(L0,Y,sum1);
            sum2 = _mm256_fmadd_pd(L1,Y,sum2);
        }
        double r0 = hsum4(sum1);
        double r1 = hsum4(sum2);
        b[i0] = b[i0] - r0;
        b[i1] = b[i1] - r1 - Lrow1[i0]*b[i0];
    }
    {
        double b3 = b[N-1];
        double b2 = b[N-2];
        double b1 = b[N-3];
        double b0 = b[N-4];
        double* __restrict__ rowN2 = &v_[(N-2)*alignN];
        double* __restrict__ rowN3 = &v_[(N-3)*alignN];
        double* __restrict__ rowN4 = &v_[(N-4)*alignN];
        const double invD33 = this->invD[N-1];
        const double invD22 = this->invD[N-2];
        const double invD11 = this->invD[N-3];
        const double invD00 = this->invD[N-4];
        double L23 = rowN2[N-1];
        double L12 = rowN3[N-2];        
        double L13 = rowN3[N-1];
        double L03 = rowN4[N-1];
        double L02 = rowN4[N-2];
        double L01 = rowN4[N-3];
        b3 = b3 * invD33;
        b2 = std::fma(-L23, b3, b2);
        b2 = b2 * invD22;
        b1 = std::fma(-L13, b3, b1);
        b1 = std::fma(-L12, b2, b1);
        b1 = b1 * invD11;
        b0 = std::fma(-L03, b3, b0);
        b0 = std::fma(-L02, b2, b0);
        b0 = std::fma(-L01, b1, b0);
        b0 = b0 * invD00;
        b[N-1] = b3;
        b[N-2] = b2;
        b[N-3] = b1;
        b[N-4] = b0;
    }

    for(int i = N-1-4; i >= 2; i=i-4)
    {
        const unsigned int i0 = i-3;
        const unsigned int i1 = i-2;
        const unsigned int i2 = i-1;
        const unsigned int i3 = i-0;
        double sum0 = b[i0];
        double sum1 = b[i1];
        double sum2 = b[i2];
        double sum3 = b[i3];
        //__m256d sumv = _mm256_loadu_pd(&b[i0]);
        __m256d result0 = _mm256_setzero_pd();
        __m256d result1 = _mm256_setzero_pd();
        __m256d result2 = _mm256_setzero_pd();
        __m256d result3 = _mm256_setzero_pd();
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ]; 
        double* __restrict__ Lrow2 = &v_[ i2 * alignN ];
        double* __restrict__ Lrow3 = &v_[ i3 * alignN ];  
        for(unsigned int j = i+1; j<this->N; j=j+4)
        {
            unsigned int j0 = j + 0;
            __m256d bv = _mm256_loadu_pd(&b[j0]);
            __m256d A00v = _mm256_loadu_pd(&Lrow0[j]);
            __m256d A10v = _mm256_loadu_pd(&Lrow1[j]);
            __m256d A20v = _mm256_loadu_pd(&Lrow2[j]);
            __m256d A30v = _mm256_loadu_pd(&Lrow3[j]);

            result0 = _mm256_fmadd_pd(A00v,bv,result0);
            result1 = _mm256_fmadd_pd(A10v,bv,result1);
            result2 = _mm256_fmadd_pd(A20v,bv,result2);
            result3 = _mm256_fmadd_pd(A30v,bv,result3);        
        }

        sum0 = sum0 - this->hsum4(result0);
        sum1 = sum1 - this->hsum4(result1);
        sum2 = sum2 - this->hsum4(result2);
        sum3 = sum3 - this->hsum4(result3);


        const double invD33 = this->invD[i3];
        const double invD22 = this->invD[i2];
        const double invD11 = this->invD[i1];
        const double invD00 = this->invD[i0];
        double L23 = Lrow2[i3];
        double L12 = Lrow1[i2];
        double L13 = Lrow1[i3];
        double L01 = Lrow0[i1];
        double L02 = Lrow0[i2];
        double L03 = Lrow0[i3];
        double b3 = b[i3];
        double b2 = b[i2];
        double b1 = b[i1];
        double b0 = b[i0];
        b3 = sum3 * invD33;
        b2 = std::fma(-L23,b3,sum2);
        b2 = b2 * invD22;
        b1 = std::fma(-L12,b2,sum1);
        b1 = std::fma(-L13,b3,b1);
        b1 = b1 * invD11;
        b0 = std::fma(-L01,b1,sum0);
        b0 = std::fma(-L02,b2,b0);
        b0 = std::fma(-L03,b3,b0);
        b0 = b0 * invD00;
        b[i3] = b3;
        b[i2] = b2;
        b[i1] = b1;
        b[i0] = b0;
    }

    {
        const unsigned int i0 = 0;
        const unsigned int i1 = 1;        
        double sum0 = b[i0];
        double sum1 = b[i1];
        __m256d result0 = _mm256_setzero_pd();
        __m256d result1 = _mm256_setzero_pd();
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ];        
        for(unsigned int j = i1+1; j<this->N; j=j+4)
        {
            unsigned int j0 = j + 0;
            __m256d bv = _mm256_loadu_pd(&b[j0]);
            __m256d A0 = _mm256_loadu_pd(&Lrow0[j]);
            __m256d A1 = _mm256_loadu_pd(&Lrow1[j]);
            result0 = _mm256_fmadd_pd(A0,bv,result0);
            result1 = _mm256_fmadd_pd(A1,bv,result1);
        }
        double r0 = hsum4(result0);
        double r1 = hsum4(result1);
        sum0 = sum0 - r0;
        sum1 = sum1 - r1;
        const double invD11 = this->invD[i1];
        const double invD00 = this->invD[i0];
        b[i1] = (sum1) * invD11;          
        b[i0] = (sum0 - Lrow0[i1]*b[i1]) * invD00;   
    }
}

void FastChemistry::LUsolver::xSolve_Vec_1
(
    double* __restrict__ b
)
{
    for(unsigned int j = 0; j<this->N;j++)
    {
        if(this->pivotIndice_[j]==j)
        {
            continue;
        }
        unsigned int jTarget = this->pivotIndice_[j];
        double temp = b[jTarget];
        b[jTarget] = b[j];
        b[j] = temp;
    }
    {
        double b0 = b[0];
        double b1 = b[1];
        double b2 = b[2];
        double b3 = b[3];
        double L10 = v_[1*alignN+0];
        double L20 = v_[2*alignN+0];
        double L21 = v_[2*alignN+1];
        double L30 = v_[3*alignN+0];
        double L31 = v_[3*alignN+1];
        double L32 = v_[3*alignN+2];
        b1 = std::fma(-L10, b0, b1);
        b2 = std::fma(-L20, b0, b2);
        b2 = std::fma(-L21, b1, b2);
        b3 = std::fma(-L30, b0, b3);
        b3 = std::fma(-L31, b1, b3);
        b3 = std::fma(-L32, b2, b3);
        b[1] = b1;
        b[2] = b2;
        b[3] = b3;
    }

    for(unsigned int i = 4; i < this->N-1; i=i+4)
    {
        unsigned int i0 = i+0;
        unsigned int i1 = i+1;
        unsigned int i2 = i+2;
        unsigned int i3 = i+3;

        __m256d sum = _mm256_loadu_pd(&b[i0]);
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ]; 
        double* __restrict__ Lrow2 = &v_[ i2 * alignN ];
        double* __restrict__ Lrow3 = &v_[ i3 * alignN ];  

        for(unsigned int j = 0; j < i; j=j+4)
        {
            __m256d L0 = _mm256_loadu_pd(&Lrow0[j]);
            __m256d L1 = _mm256_loadu_pd(&Lrow1[j]);
            __m256d L2 = _mm256_loadu_pd(&Lrow2[j]);
            __m256d L3 = _mm256_loadu_pd(&Lrow3[j]);
            __m256d Y = _mm256_loadu_pd(&b[j+0]);
            __m256d L0Y = _mm256_mul_pd(L0,Y);
            __m256d L1Y = _mm256_mul_pd(L1,Y);
            __m256d L2Y = _mm256_mul_pd(L2,Y);
            __m256d L3Y = _mm256_mul_pd(L3,Y);
            
            __m256d temp1 = _mm256_hadd_pd(L0Y,L1Y);
            __m256d temp2 = _mm256_hadd_pd(L2Y,L3Y);
            __m256d temp1_shuffle = _mm256_permute4x64_pd(temp1,_MM_SHUFFLE(3,1,2,0));
            __m256d temp2_shuffle = _mm256_permute4x64_pd(temp2,_MM_SHUFFLE(3,1,2,0));
            __m256d temp3 = _mm256_hadd_pd(temp1_shuffle,temp2_shuffle);
            __m256d temp3_shuffle = _mm256_permute4x64_pd(temp3,_MM_SHUFFLE(3,1,2,0));
            sum = _mm256_add_pd(sum,-temp3_shuffle);
        }
        _mm256_storeu_pd(&b[i0],sum);
         
        double b0 = b[i0];
        double b1 = b[i1];
        double b2 = b[i2];
        double b3 = b[i3];
        double L10 = Lrow1[i0];
        double L20 = Lrow2[i0];
        double L21 = Lrow2[i1];
        double L30 = Lrow3[i0];
        double L31 = Lrow3[i1];
        double L32 = Lrow3[i2];
        b1 = std::fma(-L10, b0, b1);
        b2 = std::fma(-L20, b0, b2);
        b2 = std::fma(-L21, b1, b2);
        b3 = std::fma(-L30, b0, b3);
        b3 = std::fma(-L31, b1, b3);
        b3 = std::fma(-L32, b2, b3);
        b[i1] = b1;
        b[i2] = b2;
        b[i3] = b3;    
    }
    {
        unsigned int i0 = this->N-1;
        int remain = i0%4;
        __m256d sum = _mm256_set1_pd(0);
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];        
        for(unsigned int j = 0; j < i0-remain; j=j+4)
        {
            __m256d L0 = _mm256_loadu_pd(&Lrow0[j]);
            __m256d Y = _mm256_loadu_pd(&b[j+0]);
            sum = _mm256_fmadd_pd(L0,Y,sum);
        }
        __m256d temp1 = _mm256_hadd_pd(sum,sum);
        __m256d temp_shuffle = _mm256_permute4x64_pd(temp1,_MM_SHUFFLE(3,1,2,0));
        __m256d temp2 = _mm256_hadd_pd(temp_shuffle,temp_shuffle);
        b[i0] = b[i0] - _mm256_cvtsd_f64(temp2);
    }
    {
        double b3 = b[N-1];
        double b2 = b[N-2];
        double b1 = b[N-3];
        double b0 = b[N-4];
        double* __restrict__ rowN2 = &v_[(N-2)*alignN];
        double* __restrict__ rowN3 = &v_[(N-3)*alignN];
        double* __restrict__ rowN4 = &v_[(N-4)*alignN];
        const double invD33 = this->invD[N-1];
        const double invD22 = this->invD[N-2];
        const double invD11 = this->invD[N-3];
        const double invD00 = this->invD[N-4];

        double L23 = rowN2[N-1];
        double L12 = rowN3[N-2];        
        double L13 = rowN3[N-1];
        double L03 = rowN4[N-1];
        double L02 = rowN4[N-2];
        double L01 = rowN4[N-3];
        b3 = b3 * invD33;
        b2 = std::fma(-L23, b3, b2);
        b2 = b2 * invD22;
        b1 = std::fma(-L13, b3, b1);
        b1 = std::fma(-L12, b2, b1);
        b1 = b1 * invD11;
        b0 = std::fma(-L03, b3, b0);
        b0 = std::fma(-L02, b2, b0);
        b0 = std::fma(-L01, b1, b0);
        b0 = b0 * invD00;
        b[N-1] = b3;
        b[N-2] = b2;
        b[N-3] = b1;
        b[N-4] = b0;    
    }


    for(int i = N-1-4; i >= 0+1; i=i-4)
    {
        const unsigned int i0 = i-3;
        const unsigned int i1 = i-2;
        const unsigned int i2 = i-1;
        const unsigned int i3 = i-0;
        double sum0 = b[i0];
        double sum1 = b[i1];
        double sum2 = b[i2];
        double sum3 = b[i3];
        __m256d result0 = _mm256_setzero_pd();
        __m256d result1 = _mm256_setzero_pd();
        __m256d result2 = _mm256_setzero_pd();
        __m256d result3 = _mm256_setzero_pd();
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ]; 
        double* __restrict__ Lrow2 = &v_[ i2 * alignN ];
        double* __restrict__ Lrow3 = &v_[ i3 * alignN ];          
        for(unsigned int j = i+1; j<this->N; j=j+4)
        {
            unsigned int j0 = j + 0;
            __m256d bv = _mm256_loadu_pd(&b[j0]);
            __m256d A00v = _mm256_loadu_pd(&Lrow0[j]);
            __m256d A10v = _mm256_loadu_pd(&Lrow1[j]);
            __m256d A20v = _mm256_loadu_pd(&Lrow2[j]);
            __m256d A30v = _mm256_loadu_pd(&Lrow3[j]);

            result0 = _mm256_fmadd_pd(A00v,bv,result0);
            result1 = _mm256_fmadd_pd(A10v,bv,result1);
            result2 = _mm256_fmadd_pd(A20v,bv,result2);
            result3 = _mm256_fmadd_pd(A30v,bv,result3);        
        }
              
        sum0 = sum0 - this->hsum4(result0);
        sum1 = sum1 - this->hsum4(result1);
        sum2 = sum2 - this->hsum4(result2);
        sum3 = sum3 - this->hsum4(result3); 
        const double invD33 = this->invD[i3];
        const double invD22 = this->invD[i2];
        const double invD11 = this->invD[i1];
        const double invD00 = this->invD[i0];

        double L23 = Lrow2[i3];
        double L12 = Lrow1[i2];
        double L13 = Lrow1[i3];
        double L01 = Lrow0[i1];
        double L02 = Lrow0[i2];
        double L03 = Lrow0[i3];

        double b3 = b[i3];
        double b2 = b[i2];
        double b1 = b[i1];
        double b0 = b[i0];

        b3 = sum3 * invD33;
        b2 = std::fma(-L23,b3,sum2);
        b2 = b2 * invD22;
        b1 = std::fma(-L12,b2,sum1);
        b1 = std::fma(-L13,b3,b1);
        b1 = b1 * invD11;
        b0 = std::fma(-L01,b1,sum0);
        b0 = std::fma(-L02,b2,b0);
        b0 = std::fma(-L03,b3,b0);
        b0 = b0 * invD00;

        b[i3] = b3;
        b[i2] = b2;
        b[i1] = b1;
        b[i0] = b0;
    }

    {
        const unsigned int i00 = 0;
        double sum00 = b[i00];
        const double invD00 = this->invD[i00];
        __m256d result0 = _mm256_setzero_pd();
        double* __restrict__ Lrow0 = &v_[ i00 * alignN ];        
        for(unsigned int j = i00+1; j<this->N; j=j+4)
        {
            unsigned int j0 = j + 0;
            __m256d A = _mm256_loadu_pd(&Lrow0[j]);
            __m256d bv = _mm256_loadu_pd(&b[j0]);
            result0 = _mm256_fmadd_pd(A,bv,result0);
 
        }
        double r0 = hsum4(result0);        
        sum00 = sum00 - r0;
        b[i00] = (sum00) * invD00;
    }
}
void FastChemistry::LUsolver::xSolve_Vec_0
(
    double* __restrict__ b
)
{

    for(unsigned int j = 0; j<this->N;j++)
    {
        if(this->pivotIndice_[j]==j)
        {
            continue;
        }
        unsigned int jTarget = this->pivotIndice_[j];
        double temp = b[jTarget];
        b[jTarget] = b[j];
        b[j] = temp;
    }
    {
        double b0 = b[0];
        double b1 = b[1];
        double b2 = b[2];
        double b3 = b[3];
        double L10 = v_[1*alignN+0];
        double L20 = v_[2*alignN+0];
        double L21 = v_[2*alignN+1];
        double L30 = v_[3*alignN+0];
        double L31 = v_[3*alignN+1];
        double L32 = v_[3*alignN+2];
        b1 = std::fma(-L10, b0, b1);
        b2 = std::fma(-L20, b0, b2);
        b2 = std::fma(-L21, b1, b2);
        b3 = std::fma(-L30, b0, b3);
        b3 = std::fma(-L31, b1, b3);
        b3 = std::fma(-L32, b2, b3);
        b[1] = b1;
        b[2] = b2;
        b[3] = b3;
    }

    for(unsigned int i = 4; i < this->N; i=i+4)
    {
        unsigned int i0 = i+0;
        unsigned int i1 = i+1;
        unsigned int i2 = i+2;
        unsigned int i3 = i+3;

        __m256d sum = _mm256_loadu_pd(&b[i0]);
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = Lrow0 + alignN; 
        double* __restrict__ Lrow2 = Lrow1 + alignN;
        double* __restrict__ Lrow3 = Lrow2 + alignN; 
        for(unsigned int j = 0; j < i; j=j+4)
        {
            __m256d L0 = _mm256_loadu_pd(&Lrow0[j]);
            __m256d L1 = _mm256_loadu_pd(&Lrow1[j]);
            __m256d L2 = _mm256_loadu_pd(&Lrow2[j]);
            __m256d L3 = _mm256_loadu_pd(&Lrow3[j]);
            __m256d Y = _mm256_loadu_pd(&b[j+0]);
            __m256d L0Y = _mm256_mul_pd(L0,Y);
            __m256d L1Y = _mm256_mul_pd(L1,Y);
            __m256d L2Y = _mm256_mul_pd(L2,Y);
            __m256d L3Y = _mm256_mul_pd(L3,Y);
            __m256d temp1 = _mm256_hadd_pd(L0Y,L1Y);
            __m256d temp2 = _mm256_hadd_pd(L2Y,L3Y);
            __m256d temp1_shuffle = _mm256_permute4x64_pd(temp1,_MM_SHUFFLE(3,1,2,0));
            __m256d temp2_shuffle = _mm256_permute4x64_pd(temp2,_MM_SHUFFLE(3,1,2,0));
            __m256d temp3 = _mm256_hadd_pd(temp1_shuffle,temp2_shuffle);
            __m256d temp3_shuffle = _mm256_permute4x64_pd(temp3,_MM_SHUFFLE(3,1,2,0));
            sum = _mm256_add_pd(sum,-temp3_shuffle);
        }
        _mm256_storeu_pd(&b[i0],sum);
        double b0 = b[i0];
        double b1 = b[i1];
        double b2 = b[i2];
        double b3 = b[i3];
        double L10 = Lrow1[i0];
        double L20 = Lrow2[i0];
        double L21 = Lrow2[i1];
        double L30 = Lrow3[i0];
        double L31 = Lrow3[i1];
        double L32 = Lrow3[i2];
        b1 = std::fma(-L10, b0, b1);
        b2 = std::fma(-L20, b0, b2);
        b2 = std::fma(-L21, b1, b2);
        b3 = std::fma(-L30, b0, b3);
        b3 = std::fma(-L31, b1, b3);
        b3 = std::fma(-L32, b2, b3);
        b[i1] = b1;
        b[i2] = b2;
        b[i3] = b3;
    }

    {
        double b3 = b[N-1];
        double b2 = b[N-2];
        double b1 = b[N-3];
        double b0 = b[N-4];
        double* __restrict__ rowN2 = &v_[(N-2)*alignN];
        double* __restrict__ rowN3 = &v_[(N-3)*alignN];
        double* __restrict__ rowN4 = &v_[(N-4)*alignN];
        const double invD33 = this->invD[N-1];
        const double invD22 = this->invD[N-2];
        const double invD11 = this->invD[N-3];
        const double invD00 = this->invD[N-4];
        double L23 = rowN2[N-1];
        double L12 = rowN3[N-2];        
        double L13 = rowN3[N-1];
        double L03 = rowN4[N-1];
        double L02 = rowN4[N-2];
        double L01 = rowN4[N-3];
        b3 = b3 * invD33;
        b2 = std::fma(-L23, b3, b2);
        b2 = b2 * invD22;
        b1 = std::fma(-L13, b3, b1);
        b1 = std::fma(-L12, b2, b1);
        b1 = b1 * invD11;
        b0 = std::fma(-L03, b3, b0);
        b0 = std::fma(-L02, b2, b0);
        b0 = std::fma(-L01, b1, b0);
        b0 = b0 * invD00;
        b[N-1] = b3;
        b[N-2] = b2;
        b[N-3] = b1;
        b[N-4] = b0;
    }
    for(int i = N-1-4; i >= 0; i=i-4)
    {
        const unsigned int i0 = i-3;
        const unsigned int i1 = i-2;
        const unsigned int i2 = i-1;
        const unsigned int i3 = i-0;
        double sum0 = b[i0];
        double sum1 = b[i1];
        double sum2 = b[i2];
        double sum3 = b[i3];
        __m256d result0 = _mm256_setzero_pd();
        __m256d result1 = _mm256_setzero_pd();
        __m256d result2 = _mm256_setzero_pd();
        __m256d result3 = _mm256_setzero_pd();
        double* __restrict__ Lrow0 = &v_[ i0 * alignN ];
        double* __restrict__ Lrow1 = &v_[ i1 * alignN ]; 
        double* __restrict__ Lrow2 = &v_[ i2 * alignN ];
        double* __restrict__ Lrow3 = &v_[ i3 * alignN ];
        for(unsigned int j = i+1; j<this->N; j=j+4)
        {
            unsigned int j0 = j + 0;
            __m256d bv = _mm256_loadu_pd(&b[j0]);
            __m256d A00v = _mm256_loadu_pd(&Lrow0[j]);
            __m256d A10v = _mm256_loadu_pd(&Lrow1[j]);
            __m256d A20v = _mm256_loadu_pd(&Lrow2[j]);
            __m256d A30v = _mm256_loadu_pd(&Lrow3[j]);
            result0 = _mm256_fmadd_pd(A00v,bv,result0);
            result1 = _mm256_fmadd_pd(A10v,bv,result1);
            result2 = _mm256_fmadd_pd(A20v,bv,result2);
            result3 = _mm256_fmadd_pd(A30v,bv,result3);
        }
        sum0 = sum0 - this->hsum4(result0);
        sum1 = sum1 - this->hsum4(result1);
        sum2 = sum2 - this->hsum4(result2);
        sum3 = sum3 - this->hsum4(result3);
        const double invD33 = this->invD[i3];
        const double invD22 = this->invD[i2];
        const double invD11 = this->invD[i1];
        const double invD00 = this->invD[i0];
        double L23 = Lrow2[i3];
        double L12 = Lrow1[i2];
        double L13 = Lrow1[i3];
        double L01 = Lrow0[i1];
        double L02 = Lrow0[i2];
        double L03 = Lrow0[i3];
        double b3 = b[i3];
        double b2 = b[i2];
        double b1 = b[i1];
        double b0 = b[i0];
        b3 = sum3*invD33;
        b2 = std::fma(-L23,b3,sum2);
        b2 = b2*invD22;
        b1 = std::fma(-L12,b2,sum1);
        b1 = std::fma(-L13,b3,b1);
        b1 = b1*invD11;
        b0 = std::fma(-L01,b1,sum0);
        b0 = std::fma(-L02,b2,b0);
        b0 = std::fma(-L03,b3,b0);
        b0 = b0*invD00;
        b[i3] = b3;
        b[i2] = b2;
        b[i1] = b1;
        b[i0] = b0;
    }
}

void FastChemistry::LUsolver::xSolve_Serial
(
    double* __restrict__ b
)
{
    for(unsigned int j = 0; j<this->N;j++)
    {
        if(this->pivotIndice_[j]==j)
        {
            continue;
        }
        unsigned int jTarget = this->pivotIndice_[j];
        double temp = b[jTarget];
        b[jTarget] = b[j];
        b[j] = temp;
    }

    for(unsigned int i = 0; i < this->N; i++)
    {
        double sum = b[i];
        for(unsigned int j = 0; j < i; j++)
        {
            sum = sum - v_[i*alignN+j]*b[j];
            
        }
        b[i] = sum;
    }

    for(int i = N-1; i >= 0; i--)
    {
        double sum = b[i];
        for(unsigned int j = i+1; j<this->N; j++)
        {
            const double xj = b[j];
            sum = sum - v_[i*alignN+j]*xj;
        }
        b[i] = sum/v_[i*alignN+i];
    }
}

void FastChemistry::LUsolver::Block4LUDecompose
(
)
{
    for(unsigned int i = 0; i < this->N;i++)
    {
        this->pivotIndice_[i] = i;
    }
    int times = (this->N-this->Remain)/4;
    if(this->Remain==0)
    {
        for(int i = 0; i < times; i ++)
        {
            for(size_t j = 0; j <N;j++)
            {
                rowPtr[j] = &v_[j*alignN];
            }
            int k0 = i*4;
            int k1 = i*4+4;
            this->LUDecompose4(k0);
            this->forwardSitituate4_0(k0,k1);
            this->permutation0(k0,k1);
            this->backSitituate4_0(k0,k1);
            //this->UpdateL22U22_Vec2_0(k0,k1);   
            int n = static_cast<int>(this->N);
            int colTail = 1 - (((n-k1-1)>>2)&1); 
            int rowTime = (n-k1)/3;
            int rowTail = (n-k1)%3;
            int colTime = (n-k1<= 4)?0:(n-k1-5)/8+1;
            //std::cout<<n-k1<<" "<<colTime<<" "<<colTail<<" "<<rowTime<<" "<<rowTail<<std::endl;
            (this->*updateTable[rowTail][colTail])(k0,k1,rowTime,colTime);
            //(this->*FuncTable[i])(k0,k1,0,0);
        }
    }
    else if(this->Remain==1)
    {
        for(int i = 0; i < times; i ++)
        {
            for(size_t j = 0; j <N;j++)
            {
                rowPtr[j] = &v_[j*alignN];
            }
            int k0 = i*4;
            int k1 = i*4+4;
            this->LUDecompose4(k0);
            this->forwardSitituate4_1(k0,k1);
            this->permutation1(k0,k1);
            this->backSitituate4_1(k0,k1);
            //this->UpdateL22U22_Vec2_1(k0,k1);
            int n = static_cast<int>(this->N);
            int colTail = 1 - (((n-k1-1)>>2)&1); 
            int rowTime = (n-k1)/3;
            int rowTail = (n-k1)%3;
            int colTime = (n-k1<= 4)?0:(n-k1-5)/8+1;
            (this->*updateTable[rowTail][colTail])(k0,k1,rowTime,colTime);
        }
        if(v_[(N-1)*alignN+N-1]==0)
        {
            v_[(N-1)*alignN+N-1] = FastChemistry::LuLimiter;
        }
        this->invD[N-1] = 1.0/v_[(N-1)*alignN+N-1];
    }
    else if(this->Remain==2)
    {
        for(int i = 0; i < times; i ++)
        {
            for(size_t j = 0; j <N;j++)
            {
                rowPtr[j] = &v_[j*alignN];
            }
            int k0 = i*4;
            int k1 = i*4+4;
            this->LUDecompose4(k0);
            this->forwardSitituate4_2(k0,k1);
            this->permutation2(k0,k1);
            this->backSitituate4_2(k0,k1);
            //this->UpdateL22U22_Vec2_2(k0,k1);   
            int n = static_cast<int>(this->N);
            int colTail = 1 - (((n-k1-1)>>2)&1); 
            int rowTime = (n-k1)/3;
            int rowTail = (n-k1)%3;
            int colTime = (n-k1<= 4)?0:(n-k1-5)/8+1;
            (this->*updateTable[rowTail][colTail])(k0,k1,rowTime,colTime);
        }
        this->LUDecompose4_2();
        this->permutation_2();     
    }
    else
    {
        for(int i = 0; i < times; i ++)
        {
            for(size_t j = 0; j <N;j++)
            {
                rowPtr[j] = &v_[j*alignN];
            }
            int k0 = i*4;
            int k1 = i*4+4;
            this->LUDecompose4(k0);
            this->forwardSitituate4_3(k0,k1);
            this->permutation3(k0,k1);
            this->backSitituate4_3(k0,k1);
            //this->UpdateL22U22_Vec2_3(k0,k1); 
            int n = static_cast<int>(this->N);
            int colTail = 1 - (((n-k1-1)>>2)&1); 
            int rowTime = (n-k1)/3;
            int rowTail = (n-k1)%3;
            int colTime = (n-k1<= 4)?0:(n-k1-5)/8+1;
            (this->*updateTable[rowTail][colTail])(k0,k1,rowTime,colTime);
        }
        this->LUDecompose4_3();
        this->permutation_3();  
    }
}

void FastChemistry::LUsolver::ReAssign(double* externalData)
{this->v_ = externalData;}

