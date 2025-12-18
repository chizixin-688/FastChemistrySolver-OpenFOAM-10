


#ifndef LUsolver_H
#define LUsolver_H
#include "Macro.H"
#include <vector>
#include <immintrin.h>
#include <string>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace FastChemistry
{
class LUsolver
{
    // Private data

        //- Matrix size
        unsigned int N;

        unsigned int alignN;

        //- Remain
        unsigned int Remain;

        //- Data pointer
        double* __restrict__ v_;

        //- Store the global index of target row.
        unsigned int* pivotIndice_;

        //- Store the row address
        std::vector<double*> rowPtr;

        //- Save the reciprocal of diagonal elements
        std::vector<double> invD;

        //- Whether the data is own
        bool owner;

        using UpdateFunc = void (LUsolver::*)(unsigned int, unsigned int, unsigned int, unsigned int);

        UpdateFunc updateTable[3][2] =
        {
            {   
                &LUsolver::UpdateL22U22_Vec_38_00,
                &LUsolver::UpdateL22U22_Vec_38_04
            },
            {   
                &LUsolver::UpdateL22U22_Vec_38_10,
                &LUsolver::UpdateL22U22_Vec_38_14
            },
            {
                &LUsolver::UpdateL22U22_Vec_38_20,
                &LUsolver::UpdateL22U22_Vec_38_24
            }
        };

public:

    // Constructor

        // Constructor from external data.
        LUsolver(double* externalData, int size);

        // Disable default constructor 
        LUsolver()=delete;  

    // Destructor

        ~LUsolver();

    // Member function

        // Access

            //- (i, j) element access.
            inline double& operator()(const unsigned int i,const unsigned int j);

            //- (i, j) element access. 
            inline double operator() (const unsigned int i,const unsigned int j)const;
 
            //- pivotIndice access.
            inline unsigned int& pivotIndice(const unsigned int i);

            //- pivotIndice access.
            inline unsigned int pivotIndice(const unsigned int i) const;

        // Print the variable

            void printMatrix(double*A,int mRows, int nCol);

            void printMatrix();

            void printPivotIndice();

        // Block LU decomposition algorithm

            // using 4*4 block.
            void Block4LUDecompose();

            // Solve Ax=b.
            void xSolve_Serial(double* __restrict__ b);
            void xSolve_Vec_0(double* __restrict__ b);
            void xSolve_Vec_1(double* __restrict__ b);
            void xSolve_Vec_2(double* __restrict__ b);
            void xSolve_Vec_3(double* __restrict__ b);
            void xSolve(double* __restrict__ b);            

            // Sub functions

                // LU decomposition of 4*4 block.
                void LUDecompose4(unsigned int k0);
                //void LUDecompose4Old(unsigned int k0);
                //void LUDecompose4();
                void LUDecompose4_2();           
                void LUDecompose4_3();        
                
                //void forwardSitituate4(unsigned int k0,unsigned int k1);
                void forwardSitituate4_0(unsigned int k0,unsigned int k1);
                void forwardSitituate4_1(unsigned int k0,unsigned int k1);
                void forwardSitituate4_2(unsigned int k0,unsigned int k1);
                void forwardSitituate4_3(unsigned int k0,unsigned int k1);

                //void backSitituate4(unsigned int k0,unsigned int k1);
                void backSitituate4_0(unsigned int k0,unsigned int k1);
                void backSitituate4_1(unsigned int k0,unsigned int k1);
                void backSitituate4_2(unsigned int k0,unsigned int k1);
                void backSitituate4_3(unsigned int k0,unsigned int k1);                

                void permutation0(unsigned int k0,unsigned int k1);
                void permutation1(unsigned int k0,unsigned int k1);
                void permutation2(unsigned int k0,unsigned int k1);
                void permutation3(unsigned int k0,unsigned int k1);
                //void permutationold(unsigned int k0,unsigned int k1);
                //void permutation();
                void permutation_2();
                void permutation_3();

                //void UpdateL22U22(unsigned int k0,unsigned int k1);
                //void UpdateL22U22_Vec(unsigned int k0,unsigned int k1);
                //void UpdateL22U22_Vec2(unsigned int k0,unsigned int k1);   
                void UpdateL22U22_Vec2_0(unsigned int k0,unsigned int k1);   
                void UpdateL22U22_Vec2_1(unsigned int k0,unsigned int k1);  
                void UpdateL22U22_Vec2_2(unsigned int k0,unsigned int k1);  
                void UpdateL22U22_Vec2_3(unsigned int k0,unsigned int k1);   

                void UpdateL22U22_Vec_38_00(unsigned int k0,unsigned int k1, unsigned int rowTimes, unsigned int colTimes);
                void UpdateL22U22_Vec_38_10(unsigned int k0,unsigned int k1, unsigned int rowTimes, unsigned int colTimes);
                void UpdateL22U22_Vec_38_20(unsigned int k0,unsigned int k1, unsigned int rowTimes, unsigned int colTimes);

                void UpdateL22U22_Vec_38_04(unsigned int k0,unsigned int k1, unsigned int rowTimes, unsigned int colTimes);
                void UpdateL22U22_Vec_38_14(unsigned int k0,unsigned int k1, unsigned int rowTimes, unsigned int colTimes);
                void UpdateL22U22_Vec_38_24(unsigned int k0,unsigned int k1, unsigned int rowTimes, unsigned int colTimes);
        //void LUDecompose_s();

        // Read square matrix from text.
        //void ReadTxt(std::string fileName);

        // Reassign a square matrix
        void ReAssign(double* externalData);

        inline double hsum4(__m256d v);// result = v0+v1+v2+v3

        inline __m256d hsum4x4(__m256d sum0, __m256d sum1,__m256d sum2, __m256d sum3) const
        {
            __m256d t0 = _mm256_hadd_pd(sum0, sum1);  // [a0+a1, b0+b1, a2+a3, b2+b3]
            __m256d t1 = _mm256_hadd_pd(sum2, sum3);  // [c0+c1, d0+d1, c2+c3, d2+d3]
            __m256d t2 = _mm256_permute4x64_pd (t0, _MM_SHUFFLE(3,1,2,0)); // [a0+a1, a2+a3, b0+b1,  b2+b3]
            __m256d t3 = _mm256_permute4x64_pd(t1, _MM_SHUFFLE(3,1,2,0)); // [c0+c1, c2+c3, d0+d1, d2+d3]
            __m256d t4 = _mm256_hadd_pd(t2, t3); // [a0+a1+a2+a3, c0+c1+c2+c3, b0+b1+b2+b3, d0+d1+d2+d3]
            __m256d t5 = _mm256_permute4x64_pd(t4, _MM_SHUFFLE(3,1,2,0));
            return t5;
        }

        inline double get_elem0(__m256d vec) const
        {
            return _mm256_cvtsd_f64(vec); 
        }

        inline double get_elem1(__m256d vec) const
        {
            __m128d low = _mm256_castpd256_pd128(vec); 
            return _mm_cvtsd_f64(_mm_unpackhi_pd(low, low));
        }

        inline double get_elem2(__m256d vec) const
        {
            __m128d high = _mm256_extractf128_pd(vec, 1); 
            return _mm_cvtsd_f64(high);
        }

        inline __m128d get_high(__m256d vec) const
        {
            return _mm256_extractf128_pd(vec, 1);
        }

        inline __m128d get_low(__m256d vec) const
        {
            return _mm256_castpd256_pd128(vec);
        }

        inline double get_elem3(__m256d vec) const
        {
            __m128d high = _mm256_extractf128_pd(vec, 1); 
            return _mm_cvtsd_f64(_mm_unpackhi_pd(high, high));
        }

        inline void transpose4x4_pd(__m256d& v0,__m256d& v1,__m256d& v2,__m256d& v3);
};

#include "LUsolverI.H"


} // End namespace Foam


//#ifdef NoRepository
//    #include "LUsolver.C" 
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
