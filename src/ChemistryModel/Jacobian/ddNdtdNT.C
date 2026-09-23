/*---------------------------------------------------------------------------*\
  Description
      Jacobian-matrix assembly in molar-amount (N) space.

      Reaction rates are obtained as functions of molar concentrations (c)
      and temperature (T).  The routines below convert the concentration
      based Jacobian d(dN/dt)/d(c,T), supplied by the reaction mechanism,
      into the molar-amount based Jacobian consumed by the ODE solver:

          ddNTdtdNT[i*nCols+j]  = d(dNi/dt)/dNj          species block
          ddNTdtdNT[i*nCols+Ns]  = d(dNi/dt)/dT          last column (T)
          ddNTdtdNT[Ns*nCols+j]  = d(dT/dt)/dNj          last row (T)
          ddNTdtdNT[Ns*nCols+Ns] = d(dT/dt)/dT

      and dNTdt is converted from d(c)/dt to dN/dt.  The tail is zero.

      Scalar reference implementations (ddNdtdNT, ddTdtdNT, ddNTdtdNT) are
      provided first, followed by the AVX2 SIMD kernels specialised for the
      species-count remainder Ns % 4 == 0..3 (ddNTdtdNT0..3).  The
      '...WithMergeI' variants merge the Jacobian with the identity shift
      required by the implicit ODE matrix (A = J - I*scale).

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "FastChemistryModel.H"
#include "vec_math_avx2.H"

//---------------------------------
// 2. SIMD / AVX2 headers
//---------------------------------
#include <immintrin.h>

namespace Foam
{

//=============================================================================//
//  Scalar reference implementations (generic Ns)
//=============================================================================//

// ddNdtdNT: species block  d(dN_i/dt)/dN_j  and the temperature column
// d(dN_i/dt)/dT (column Ns); also converts dNTdt from d(cV)/dt to dN/dt.

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNdtdNT
(
    const double* __restrict__ Phi,
    const double* __restrict__ c,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    int Ns,
    int nCols
) const noexcept
{
    for (int i=0; i<Ns; i++)
    {
        double sum = 0;
        for (int k=0; k<Ns; k++)
        {
            const double ddNidtByVdck = ddNTdtdNT[i*nCols+k];
            sum += ddNidtByVdck*c[k];
        }

        for (int j=0; j<Ns; j++)
        {
            const double omegai = dNTdt[i];
            ddNTdtdNT[i*nCols+j] = (omegai-sum)*V*invNtot + ddNTdtdNT[i*nCols+j];
        }

        // Convert dcdt to dNidt
        dNTdt[i] = dNTdt[i]*V;

        double ddNidtByVdT = ddNTdtdNT[i*nCols+Ns] - sum*invT;
        double& ddNidtdT = ddNTdtdNT[i*nCols+Ns];
        ddNidtdT = V*ddNidtByVdT + dNTdt[i]*invT;
    }
}

// ddTdtdNT: temperature row  d(dT/dt)/dN_j  and the diagonal entry
// d(dT/dt)/dT (row Ns), based on the species rates in dNTdt/ddNTdtdNT.

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddTdtdNT
(
    double* __restrict__ dNTdt,
    const double* __restrict__ NT,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ ddNTdtdNT,
    double invcptot,
    double dcptotdT,
    int Ns,
    int cols
) const noexcept
{
    // dT/dt
    double& dTdt = dNTdt[Ns];
    dTdt = 0;
    for (label i=0; i<Ns; i++)
    {
        dTdt -= dNTdt[i]*ha[i];
    }
    dTdt = dTdt*invcptot;

    // d(dTdt)/dY
    for (label i=0; i<Ns; i++)
    {
        double& ddTdtdni = ddNTdtdNT[Ns*cols+i];
        ddTdtdni = 0;
        for (label j=0; j<Ns; j++)
        {
            const double ddnjdtdni = ddNTdtdNT[j*cols+i];
            ddTdtdni -= ddnjdtdni*ha[j];
        }
        ddTdtdni -= cp[i]*dTdt;
        ddTdtdni = ddTdtdni*invcptot;
    }

    // d(dTdt)/dT
    double& ddTdtdT = ddNTdtdNT[Ns*cols+Ns];
    ddTdtdT = 0;
    for (label i=0; i<Ns; i++)
    {
        const double ddNidtdT = ddNTdtdNT[i*cols+Ns];
        ddTdtdT -= dNTdt[i]*cp[i] + ddNidtdT*ha[i];
    }
    ddTdtdT -= dTdt*dcptotdT;
    ddTdtdT = ddTdtdT*invcptot;
}

// ddNTdtdNT: fused scalar implementation that assembles the complete
// Jacobian (species block + T column + T row) in a single pass.

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNTdtdNT
(
    const double* __restrict__ NT,
    const double* __restrict__ c,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    double invcptot,
    double dcptotdT,
    int Ns,
    int nCols
) const noexcept
{
    // dT/dt
    double& dTdt = dNTdt[Ns];
    dTdt = 0;
    for (int i=0; i<Ns; i++)
    {
        dTdt -= dNTdt[i]*ha[i]*V;
        ddNTdtdNT[Ns*nCols+i] = 0;
    }
    dTdt = dTdt*invcptot;

    for (int i=0; i<Ns; i++)
    {
        double sum = 0;
        for (int k=0; k<Ns; k++)
        {
            const double ddNidtByVdck = ddNTdtdNT[i*nCols+k];
            sum += ddNidtByVdck*c[k];
        }

        for (int j=0; j<Ns; j++)
        {
            const double omegai = dNTdt[i];
            const double ddNidtdNj = (omegai-sum)*V*invNtot + ddNTdtdNT[i*nCols+j];
            ddNTdtdNT[i*nCols+j] =  ddNidtdNj;
            ddNTdtdNT[Ns*nCols+j] -= ddNidtdNj*ha[i];
        }

        // Convert dcdt to dNidt
        dNTdt[i] = dNTdt[i]*V;

        double ddNidtByVdT = ddNTdtdNT[i*nCols+Ns] - sum*invT;

        const double ddNidtdT = V*ddNidtByVdT + dNTdt[i]*invT;
        ddNTdtdNT[i*nCols+Ns] = ddNidtdT;
        ddNTdtdNT[Ns*nCols+Ns] -= dNTdt[i]*cp[i] + ddNidtdT*ha[i];
    }

    for (label i=0; i<Ns; i++)
    {
        ddNTdtdNT[Ns*nCols+i] -= cp[i]*dTdt;
        ddNTdtdNT[Ns*nCols+i] *= invcptot;
    }
    ddNTdtdNT[Ns*nCols+Ns] -= dTdt*dcptotdT;
    ddNTdtdNT[Ns*nCols+Ns] = ddNTdtdNT[Ns*nCols+Ns]*invcptot;
}

//=============================================================================//
//  AVX2 SIMD implementations specialised for Ns % 4 == 0..3
//
//      ddNTdtdNT0..3            : full Jacobian assembly
//      ddNTdtdNT0..3WithMergeI  : fused Jacobian + I*scale assembly
//=============================================================================//

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNTdtdNT0
(
    const double* __restrict__ NT,
    const double* __restrict__ c,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    double invcptot,
    double dcptotdT,
    int Ns,
    int nCols
) const noexcept
{
    for (int i=0; i<Ns; i++)
    {
        ddNTdtdNT[Ns*nCols+i] = 0;
    }
    // dT/dt
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<Ns; i=i+4)
    {
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d(&dNTdt[i+0]), load256d(&ha[i+0])));
    }
    dNTdt[Ns] = hsum4(dTdtv)*V*invcptot;
    __m256d VByNtotv = _mm256_set1_pd(V*invNtot);
    __m256d Vv = _mm256_set1_pd(V);

    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (int i=0; i<Ns; i=i+4)
    {

        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        for (int k=0; k<Ns; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi3v);
        __m256d dNi03dt = load256d(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d result3v = _mm256_permute4x64_pd(result,0b11111111);
        __m256d hai03v = load256d(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        __m256d hai3v = _mm256_set1_pd(ha[i+3]);
        for (int j=0; j<Ns; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],Jaci2j03v);
            store256d(&ddNTdtdNT[(i+3)*nCols+j+0],Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],ddNTdtdNT[(i+3)*nCols+Ns]);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = get2(ddNidtdTv);
        ddNTdtdNT[(i+3)*nCols+Ns] = get3(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d(&cp[i+0]),ddTdtdTv);
    }
    const double dTdt = dNTdt[Ns];
    dTdtv = _mm256_set1_pd(dNTdt[Ns]);
    __m256d invcptotv = _mm256_set1_pd(invcptot);
    for (label i=0; i<Ns; i=i+4)
    {
        __m256d ddTdtdNi0v = fnmadd256d(load256d(&cp[i+0]),dTdtv,load256d(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d(&ddNTdtdNT[Ns*nCols+i],ddTdtdNi0v);
    }
    double ddTdtdT = hsum4(ddTdtdTv);
    ddTdtdT = ddTdtdT - dTdt*dcptotdT;
    ddTdtdT = ddTdtdT*invcptot;
    ddNTdtdNT[Ns*nCols+Ns] = ddTdtdT;
}

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNTdtdNT0WithMergeI
(
    const double* __restrict__ NT,
    const double* __restrict__ c,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    double invcptot,
    double dcptotdT,
    double scale,
    int Ns,
    int nCols
) const noexcept
{
    for (int i=0; i<Ns; i++)
    {
        ddNTdtdNT[Ns*nCols+i] = 0;
    }
    // dT/dt
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<Ns; i=i+4)
    {
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d(&dNTdt[i+0]), load256d(&ha[i+0])));
    }
    dNTdt[Ns] = hsum4(dTdtv)*V*invcptot;
    __m256d VByNtotv = _mm256_set1_pd(V*invNtot);
    __m256d Vv = _mm256_set1_pd(V);

    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (int i=0; i<Ns; i=i+4)
    {

        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        for (int k=0; k<Ns; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi3v);
        __m256d dNi03dt = load256d(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d result3v = _mm256_permute4x64_pd(result,0b11111111);
        __m256d hai03v = load256d(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        __m256d hai3v = _mm256_set1_pd(ha[i+3]);
        for (int j=0; j<Ns; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],-Jaci2j03v);
            store256d(&ddNTdtdNT[(i+3)*nCols+j+0],-Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        ddNTdtdNT[(i+0)*nCols+i+0] += scale;
        ddNTdtdNT[(i+1)*nCols+i+1] += scale;
        ddNTdtdNT[(i+2)*nCols+i+2] += scale;
        ddNTdtdNT[(i+3)*nCols+i+3] += scale;
        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],ddNTdtdNT[(i+3)*nCols+Ns]);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = -get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = -get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = -get2(ddNidtdTv);
        ddNTdtdNT[(i+3)*nCols+Ns] = -get3(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d(&cp[i+0]),ddTdtdTv);
    }
    const double dTdt = dNTdt[Ns];
    dTdtv = _mm256_set1_pd(dNTdt[Ns]);
    __m256d invcptotv = _mm256_set1_pd(invcptot);
    for (label i=0; i<Ns; i=i+4)
    {
        __m256d ddTdtdNi0v = fnmadd256d(load256d(&cp[i+0]),dTdtv,load256d(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d(&ddNTdtdNT[Ns*nCols+i],-ddTdtdNi0v);
    }
    double ddTdtdT = hsum4(ddTdtdTv);
    ddTdtdT = ddTdtdT - dTdt*dcptotdT;
    ddTdtdT = ddTdtdT*invcptot;
    ddNTdtdNT[Ns*nCols+Ns] = -ddTdtdT + scale;
}
template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNTdtdNT1
(
    const double* __restrict__ NT,
    const double* __restrict__ c,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    double invcptot,
    double dcptotdT,
    int Ns,
    int nCols
) const noexcept
{
    for (int i=0; i<Ns; i++)
    {
        ddNTdtdNT[Ns*nCols+i] = 0;
    }
    // dT/dt
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<Ns-1; i=i+4)
    {
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d(&dNTdt[i+0]), load256d(&ha[i+0])));
    }
    {
        int i=Ns-1;
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d0(&dNTdt[i+0]), load256d0(&ha[i+0])));
    }
    dNTdt[Ns] = hsum4(dTdtv)*V*invcptot;
    __m256d VByNtotv = _mm256_set1_pd(V*invNtot);
    __m256d Vv = _mm256_set1_pd(V);

    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (int i=0; i<Ns-1; i=i+4)
    {

        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        for (int k=0; k<Ns-1; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }
        {
            int k=Ns-1;
            __m256d c03v = load256d0(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi3v);
        __m256d dNi03dt = load256d(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d result3v = _mm256_permute4x64_pd(result,0b11111111);
        __m256d hai03v = load256d(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        __m256d hai3v = _mm256_set1_pd(ha[i+3]);
        for (int j=0; j<Ns-1; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],Jaci2j03v);
            store256d(&ddNTdtdNT[(i+3)*nCols+j+0],Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-1;
            __m256d Jaci0j03v = add256d(result0v,load256d0(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d0(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d0(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d0(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d0(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d0(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);
            store256d0(&ddNTdtdNT[(i+2)*nCols+j+0],Jaci2j03v);
            store256d0(&ddNTdtdNT[(i+3)*nCols+j+0],Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d0(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d0(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],ddNTdtdNT[(i+3)*nCols+Ns]);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = get2(ddNidtdTv);
        ddNTdtdNT[(i+3)*nCols+Ns] = get3(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d(&cp[i+0]),ddTdtdTv);
    }
    {
        int i = Ns-1;
        __m256d sumi0v = _mm256_setzero_pd();
        for (int k=0; k<Ns-1; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
        }
        {
            int k=Ns-1;
            __m256d c03v = load256d0(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
        }

        __m256d sumi03=  _mm256_setr_pd(hsum4(sumi0v),0,0,0);
        __m256d dNi03dt = load256d0(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d hai03v = load256d0(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        for (int j=0; j<Ns-1; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-1;
            __m256d Jaci0j03v = add256d(result0v,load256d0(&ddNTdtdNT[(i+0)*nCols+j+0]));

            store256d0(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d0(&ddNTdtdNT[Ns*nCols+j]));
            store256d0(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d0(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d0(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],0,0,0);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = get0(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d0(&cp[i+0]),ddTdtdTv);
    }
    const double dTdt = dNTdt[Ns];
    dTdtv = _mm256_set1_pd(dNTdt[Ns]);
    __m256d invcptotv = _mm256_set1_pd(invcptot);
    for (int i=0; i<Ns-1; i=i+4)
    {
        __m256d ddTdtdNi0v = fnmadd256d(load256d(&cp[i+0]),dTdtv,load256d(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d(&ddNTdtdNT[Ns*nCols+i],ddTdtdNi0v);
    }
    {
        int i=Ns-1;

        __m256d ddTdtdNi0v = fnmadd256d(load256d0(&cp[i+0]),dTdtv,load256d0(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d0(&ddNTdtdNT[Ns*nCols+i],ddTdtdNi0v);
    }
    double ddTdtdT = hsum4(ddTdtdTv);
    ddTdtdT = ddTdtdT - dTdt*dcptotdT;
    ddTdtdT = ddTdtdT*invcptot;
    ddNTdtdNT[Ns*nCols+Ns] = ddTdtdT;
}

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNTdtdNT1WithMergeI
(
    const double* __restrict__ NT,
    const double* __restrict__ c,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    double invcptot,
    double dcptotdT,
    double scale,
    int Ns,
    int nCols
) const noexcept
{
    for (int i=0; i<Ns; i++)
    {
        ddNTdtdNT[Ns*nCols+i] = 0;
    }
    // dT/dt
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<Ns-1; i=i+4)
    {
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d(&dNTdt[i+0]), load256d(&ha[i+0])));
    }
    {
        int i=Ns-1;
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d0(&dNTdt[i+0]), load256d0(&ha[i+0])));
    }
    dNTdt[Ns] = hsum4(dTdtv)*V*invcptot;
    __m256d VByNtotv = _mm256_set1_pd(V*invNtot);
    __m256d Vv = _mm256_set1_pd(V);

    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (int i=0; i<Ns-1; i=i+4)
    {

        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        for (int k=0; k<Ns-1; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }
        {
            int k=Ns-1;
            __m256d c03v = load256d0(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi3v);
        __m256d dNi03dt = load256d(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d result3v = _mm256_permute4x64_pd(result,0b11111111);
        __m256d hai03v = load256d(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        __m256d hai3v = _mm256_set1_pd(ha[i+3]);
        for (int j=0; j<Ns-1; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],-Jaci2j03v);
            store256d(&ddNTdtdNT[(i+3)*nCols+j+0],-Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-1;
            __m256d Jaci0j03v = add256d(result0v,load256d0(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d0(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d0(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d0(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d0(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d0(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);
            store256d0(&ddNTdtdNT[(i+2)*nCols+j+0],-Jaci2j03v);
            store256d0(&ddNTdtdNT[(i+3)*nCols+j+0],-Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d0(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d0(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        ddNTdtdNT[(i+0)*nCols+i+0] += scale;
        ddNTdtdNT[(i+1)*nCols+i+1] += scale;
        ddNTdtdNT[(i+2)*nCols+i+2] += scale;
        ddNTdtdNT[(i+3)*nCols+i+3] += scale;
        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],ddNTdtdNT[(i+3)*nCols+Ns]);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = -get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = -get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = -get2(ddNidtdTv);
        ddNTdtdNT[(i+3)*nCols+Ns] = -get3(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d(&cp[i+0]),ddTdtdTv);
    }
    {
        int i = Ns-1;
        __m256d sumi0v = _mm256_setzero_pd();
        for (int k=0; k<Ns-1; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
        }
        {
            int k=Ns-1;
            __m256d c03v = load256d0(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d0(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
        }

        __m256d sumi03=  _mm256_setr_pd(hsum4(sumi0v),0,0,0);
        __m256d dNi03dt = load256d0(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d hai03v = load256d0(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        for (int j=0; j<Ns-1; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-1;
            __m256d Jaci0j03v = add256d(result0v,load256d0(&ddNTdtdNT[(i+0)*nCols+j+0]));

            store256d0(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d0(&ddNTdtdNT[Ns*nCols+j]));
            store256d0(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        ddNTdtdNT[(i+0)*nCols+i+0] += scale;
        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d0(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d0(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],0,0,0);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = -get0(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d0(&cp[i+0]),ddTdtdTv);
    }
    const double dTdt = dNTdt[Ns];
    dTdtv = _mm256_set1_pd(dNTdt[Ns]);
    __m256d invcptotv = _mm256_set1_pd(invcptot);
    for (int i=0; i<Ns-1; i=i+4)
    {
        __m256d ddTdtdNi0v = fnmadd256d(load256d(&cp[i+0]),dTdtv,load256d(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d(&ddNTdtdNT[Ns*nCols+i],-ddTdtdNi0v);
    }
    {
        int i=Ns-1;

        __m256d ddTdtdNi0v = fnmadd256d(load256d0(&cp[i+0]),dTdtv,load256d0(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d0(&ddNTdtdNT[Ns*nCols+i],-ddTdtdNi0v);
    }
    double ddTdtdT = hsum4(ddTdtdTv);
    ddTdtdT = ddTdtdT - dTdt*dcptotdT;
    ddTdtdT = ddTdtdT*invcptot;
    ddNTdtdNT[Ns*nCols+Ns] = -ddTdtdT + scale;
}

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNTdtdNT2
(
    const double* __restrict__ NT,
    const double* __restrict__ c,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    double invcptot,
    double dcptotdT,
    int Ns,
    int nCols
) const noexcept
{
    for (int i=0; i<Ns; i++)
    {
        ddNTdtdNT[Ns*nCols+i] = 0;
    }
    // dT/dt
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<Ns-2; i=i+4)
    {
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d(&dNTdt[i+0]), load256d(&ha[i+0])));
    }
    {
        int i=Ns-2;
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d01(&dNTdt[i+0]), load256d01(&ha[i+0])));
    }
    dNTdt[Ns] = hsum4(dTdtv)*V*invcptot;
    __m256d VByNtotv = _mm256_set1_pd(V*invNtot);
    __m256d Vv = _mm256_set1_pd(V);

    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (int i=0; i<Ns-2; i=i+4)
    {
        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        for (int k=0; k<Ns-2; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }
        {
            int k=Ns-2;
            __m256d c03v = load256d01(&c[k+0]);
            sumi0v = fmadd256d(load256d01(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = fmadd256d(load256d01(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = fmadd256d(load256d01(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = fmadd256d(load256d01(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi3v);
        __m256d dNi03dt = load256d(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d result3v = _mm256_permute4x64_pd(result,0b11111111);
        __m256d hai03v = load256d(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        __m256d hai3v = _mm256_set1_pd(ha[i+3]);

        for (int j=0; j<Ns-2; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],Jaci2j03v);
            store256d(&ddNTdtdNT[(i+3)*nCols+j+0],Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-2;
            __m256d Jaci0j03v = add256d(result0v,load256d01(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d01(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d01(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d01(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d01(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d01(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);
            store256d01(&ddNTdtdNT[(i+2)*nCols+j+0],Jaci2j03v);
            store256d01(&ddNTdtdNT[(i+3)*nCols+j+0],Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d01(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d01(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],ddNTdtdNT[(i+3)*nCols+Ns]);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = get2(ddNidtdTv);
        ddNTdtdNT[(i+3)*nCols+Ns] = get3(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d(&cp[i+0]),ddTdtdTv);
    }

    {
        // tail: last two rows Ns-2, Ns-1
        int i = Ns-2;
        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        for (int k=0; k<Ns-2; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
        }
        {
            int k=Ns-2;
            __m256d c03v = load256d01(&c[k+0]);
            sumi0v = fmadd256d(load256d01(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = fmadd256d(load256d01(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi0v,sumi1v);
        __m256d dNi03dt = load256d01(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d hai03v = load256d01(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        for (int j=0; j<Ns-2; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-2;
            __m256d Jaci0j03v = add256d(result0v,load256d01(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d01(&ddNTdtdNT[(i+1)*nCols+j+0]));

            store256d01(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d01(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d01(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            store256d01(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d01(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d01(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],0,0);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = get1(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d01(&cp[i+0]),ddTdtdTv);
    }

    const double dTdt = dNTdt[Ns];
    dTdtv = _mm256_set1_pd(dNTdt[Ns]);
    __m256d invcptotv = _mm256_set1_pd(invcptot);
    for (int i=0; i<Ns-2; i=i+4)
    {
        __m256d ddTdtdNi0v = fnmadd256d(load256d(&cp[i+0]),dTdtv,load256d(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d(&ddNTdtdNT[Ns*nCols+i],ddTdtdNi0v);
    }
    {
        int i=Ns-2;
        __m256d ddTdtdNi0v = fnmadd256d(load256d01(&cp[i+0]),dTdtv,load256d01(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d01(&ddNTdtdNT[Ns*nCols+i],ddTdtdNi0v);
    }
    double ddTdtdT = hsum4(ddTdtdTv);
    ddTdtdT = ddTdtdT - dTdt*dcptotdT;
    ddTdtdT = ddTdtdT*invcptot;
    ddNTdtdNT[Ns*nCols+Ns] = ddTdtdT;
}

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNTdtdNT2WithMergeI
(
    const double* __restrict__ NT,
    const double* __restrict__ c,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    double invcptot,
    double dcptotdT,
    double scale,
    int Ns,
    int nCols
) const noexcept
{
    for (int i=0; i<Ns; i++)
    {
        ddNTdtdNT[Ns*nCols+i] = 0;
    }
    // dT/dt
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<Ns-2; i=i+4)
    {
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d(&dNTdt[i+0]), load256d(&ha[i+0])));
    }
    {
        int i=Ns-2;
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d01(&dNTdt[i+0]), load256d01(&ha[i+0])));
    }
    dNTdt[Ns] = hsum4(dTdtv)*V*invcptot;
    __m256d VByNtotv = _mm256_set1_pd(V*invNtot);
    __m256d Vv = _mm256_set1_pd(V);

    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (int i=0; i<Ns-2; i=i+4)
    {
        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        for (int k=0; k<Ns-2; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }
        {
            int k=Ns-2;
            __m256d c03v = load256d01(&c[k+0]);
            sumi0v = fmadd256d(load256d01(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = fmadd256d(load256d01(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = fmadd256d(load256d01(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = fmadd256d(load256d01(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi3v);
        __m256d dNi03dt = load256d(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d result3v = _mm256_permute4x64_pd(result,0b11111111);
        __m256d hai03v = load256d(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        __m256d hai3v = _mm256_set1_pd(ha[i+3]);

        for (int j=0; j<Ns-2; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],-Jaci2j03v);
            store256d(&ddNTdtdNT[(i+3)*nCols+j+0],-Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-2;
            __m256d Jaci0j03v = add256d(result0v,load256d01(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d01(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d01(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d01(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d01(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d01(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);
            store256d01(&ddNTdtdNT[(i+2)*nCols+j+0],-Jaci2j03v);
            store256d01(&ddNTdtdNT[(i+3)*nCols+j+0],-Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d01(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d01(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        ddNTdtdNT[(i+0)*nCols+i+0] +=scale;
        ddNTdtdNT[(i+1)*nCols+i+1] +=scale;
        ddNTdtdNT[(i+2)*nCols+i+2] +=scale;
        ddNTdtdNT[(i+3)*nCols+i+3] +=scale;
        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],ddNTdtdNT[(i+3)*nCols+Ns]);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = -get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = -get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = -get2(ddNidtdTv);
        ddNTdtdNT[(i+3)*nCols+Ns] = -get3(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d(&cp[i+0]),ddTdtdTv);
    }

    {
        // tail: last two rows Ns-2, Ns-1
        int i = Ns-2;
        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        for (int k=0; k<Ns-2; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
        }
        {
            int k=Ns-2;
            __m256d c03v = load256d01(&c[k+0]);
            sumi0v = fmadd256d(load256d01(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = fmadd256d(load256d01(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi0v,sumi1v);
        __m256d dNi03dt = load256d01(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d hai03v = load256d01(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        for (int j=0; j<Ns-2; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-2;
            __m256d Jaci0j03v = add256d(result0v,load256d01(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d01(&ddNTdtdNT[(i+1)*nCols+j+0]));

            store256d01(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d01(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d01(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            store256d01(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        ddNTdtdNT[(i+0)*nCols+i+0] += scale;
        ddNTdtdNT[(i+1)*nCols+i+1] += scale;
        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d01(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d01(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],0,0);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = -get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = -get1(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d01(&cp[i+0]),ddTdtdTv);
    }

    const double dTdt = dNTdt[Ns];
    dTdtv = _mm256_set1_pd(dNTdt[Ns]);
    __m256d invcptotv = _mm256_set1_pd(invcptot);
    for (int i=0; i<Ns-2; i=i+4)
    {
        __m256d ddTdtdNi0v = fnmadd256d(load256d(&cp[i+0]),dTdtv,load256d(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d(&ddNTdtdNT[Ns*nCols+i],-ddTdtdNi0v);
    }
    {
        int i=Ns-2;
        __m256d ddTdtdNi0v = fnmadd256d(load256d01(&cp[i+0]),dTdtv,load256d01(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d01(&ddNTdtdNT[Ns*nCols+i],-ddTdtdNi0v);
    }
    double ddTdtdT = hsum4(ddTdtdTv);
    ddTdtdT = ddTdtdT - dTdt*dcptotdT;
    ddTdtdT = ddTdtdT*invcptot;
    ddNTdtdNT[Ns*nCols+Ns] = -ddTdtdT + scale;
}

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNTdtdNT3
(
    const double* __restrict__ NT,
    const double* __restrict__ c,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    double invcptot,
    double dcptotdT,
    int Ns,
    int nCols
) const noexcept
{
    for (int i=0; i<Ns; i++)
    {
        ddNTdtdNT[Ns*nCols+i] = 0;
    }
    // dT/dt
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<Ns-3; i=i+4)
    {
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d(&dNTdt[i+0]), load256d(&ha[i+0])));
    }
    {
        int i=Ns-3;
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d012(&dNTdt[i+0]), load256d012(&ha[i+0])));
    }
    dNTdt[Ns] = hsum4(dTdtv)*V*invcptot;
    __m256d VByNtotv = _mm256_set1_pd(V*invNtot);
    __m256d Vv = _mm256_set1_pd(V);

    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (int i=0; i<Ns-3; i=i+4)
    {
        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        for (int k=0; k<Ns-3; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }
        {
            int k=Ns-3;
            __m256d c03v = load256d012(&c[k+0]);
            sumi0v = fmadd256d(load256d012(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = fmadd256d(load256d012(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = fmadd256d(load256d012(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = fmadd256d(load256d012(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi3v);
        __m256d dNi03dt = load256d(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d result3v = _mm256_permute4x64_pd(result,0b11111111);
        __m256d hai03v = load256d(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        __m256d hai3v = _mm256_set1_pd(ha[i+3]);

        for (int j=0; j<Ns-3; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],Jaci2j03v);
            store256d(&ddNTdtdNT[(i+3)*nCols+j+0],Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-3;
            __m256d Jaci0j03v = add256d(result0v,load256d012(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d012(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d012(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d012(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d012(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d012(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);
            store256d012(&ddNTdtdNT[(i+2)*nCols+j+0],Jaci2j03v);
            store256d012(&ddNTdtdNT[(i+3)*nCols+j+0],Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d012(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d012(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],ddNTdtdNT[(i+3)*nCols+Ns]);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = get2(ddNidtdTv);
        ddNTdtdNT[(i+3)*nCols+Ns] = get3(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d(&cp[i+0]),ddTdtdTv);
    }

    {
        // tail: last three rows Ns-3, Ns-2, Ns-1
        int i = Ns-3;
        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        for (int k=0; k<Ns-3; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
        }
        {
            int k=Ns-3;
            __m256d c03v = load256d012(&c[k+0]);
            sumi0v = fmadd256d(load256d012(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = fmadd256d(load256d012(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = fmadd256d(load256d012(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi2v);
        __m256d dNi03dt = load256d012(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d hai03v = load256d012(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        for (int j=0; j<Ns-3; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],Jaci2j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-3;
            __m256d Jaci0j03v = add256d(result0v,load256d012(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d012(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d012(&ddNTdtdNT[(i+2)*nCols+j+0]));

            store256d012(&ddNTdtdNT[(i+0)*nCols+j+0],Jaci0j03v);
            store256d012(&ddNTdtdNT[(i+1)*nCols+j+0],Jaci1j03v);
            store256d012(&ddNTdtdNT[(i+2)*nCols+j+0],Jaci2j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d012(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            store256d012(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d012(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d012(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],0);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = get2(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d012(&cp[i+0]),ddTdtdTv);
    }

    const double dTdt = dNTdt[Ns];
    dTdtv = _mm256_set1_pd(dNTdt[Ns]);
    __m256d invcptotv = _mm256_set1_pd(invcptot);
    for (int i=0; i<Ns-3; i=i+4)
    {
        __m256d ddTdtdNi0v = fnmadd256d(load256d(&cp[i+0]),dTdtv,load256d(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d(&ddNTdtdNT[Ns*nCols+i],ddTdtdNi0v);
    }
    {
        int i=Ns-3;
        __m256d ddTdtdNi0v = fnmadd256d(load256d012(&cp[i+0]),dTdtv,load256d012(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d012(&ddNTdtdNT[Ns*nCols+i],ddTdtdNi0v);
    }
    double ddTdtdT = hsum4(ddTdtdTv);
    ddTdtdT = ddTdtdT - dTdt*dcptotdT;
    ddTdtdT = ddTdtdT*invcptot;
    ddNTdtdNT[Ns*nCols+Ns] = ddTdtdT;
}

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddNTdtdNT3WithMergeI
(
    const double* __restrict__ NT,
    const double* __restrict__ c,
    const double* __restrict__ ha,
    const double* __restrict__ cp,
    double* __restrict__ dNTdt,
    double* __restrict__ ddNTdtdNT,
    double invNtot,
    double invT,
    double V,
    double invcptot,
    double dcptotdT,
    double scale,
    int Ns,
    int nCols
) const noexcept
{
    for (int i=0; i<Ns; i++)
    {
        ddNTdtdNT[Ns*nCols+i] = 0;
    }
    // dT/dt
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<Ns-3; i=i+4)
    {
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d(&dNTdt[i+0]), load256d(&ha[i+0])));
    }
    {
        int i=Ns-3;
        dTdtv = _mm256_sub_pd(dTdtv, _mm256_mul_pd(load256d012(&dNTdt[i+0]), load256d012(&ha[i+0])));
    }
    dNTdt[Ns] = hsum4(dTdtv)*V*invcptot;
    __m256d VByNtotv = _mm256_set1_pd(V*invNtot);
    __m256d Vv = _mm256_set1_pd(V);

    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (int i=0; i<Ns-3; i=i+4)
    {
        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        for (int k=0; k<Ns-3; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }
        {
            int k=Ns-3;
            __m256d c03v = load256d012(&c[k+0]);
            sumi0v = fmadd256d(load256d012(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = fmadd256d(load256d012(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = fmadd256d(load256d012(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
            sumi3v = fmadd256d(load256d012(&ddNTdtdNT[(i+3)*nCols+k+0]),c03v,sumi3v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi3v);
        __m256d dNi03dt = load256d(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d result3v = _mm256_permute4x64_pd(result,0b11111111);
        __m256d hai03v = load256d(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        __m256d hai3v = _mm256_set1_pd(ha[i+3]);

        for (int j=0; j<Ns-3; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],-Jaci2j03v);
            store256d(&ddNTdtdNT[(i+3)*nCols+j+0],-Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-3;
            __m256d Jaci0j03v = add256d(result0v,load256d012(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d012(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d012(&ddNTdtdNT[(i+2)*nCols+j+0]));
            __m256d Jaci3j03v = add256d(result3v,load256d012(&ddNTdtdNT[(i+3)*nCols+j+0]));

            store256d012(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d012(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);
            store256d012(&ddNTdtdNT[(i+2)*nCols+j+0],-Jaci2j03v);
            store256d012(&ddNTdtdNT[(i+3)*nCols+j+0],-Jaci3j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d012(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai3v,Jaci3j03v,ddTdtdNj);
            store256d012(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        ddNTdtdNT[(i+0)*nCols+i+0] += scale;
        ddNTdtdNT[(i+1)*nCols+i+1] += scale;
        ddNTdtdNT[(i+2)*nCols+i+2] += scale;
        ddNTdtdNT[(i+3)*nCols+i+3] += scale;
        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],ddNTdtdNT[(i+3)*nCols+Ns]);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = -get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = -get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = -get2(ddNidtdTv);
        ddNTdtdNT[(i+3)*nCols+Ns] = -get3(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d(&cp[i+0]),ddTdtdTv);
    }

    {
        // tail: last three rows Ns-3, Ns-2, Ns-1
        int i = Ns-3;
        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        for (int k=0; k<Ns-3; k=k+4)
        {
            __m256d c03v = load256d(&c[k+0]);
            sumi0v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = _mm256_fmadd_pd(load256d(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
        }
        {
            int k=Ns-3;
            __m256d c03v = load256d012(&c[k+0]);
            sumi0v = fmadd256d(load256d012(&ddNTdtdNT[(i+0)*nCols+k+0]),c03v,sumi0v);
            sumi1v = fmadd256d(load256d012(&ddNTdtdNT[(i+1)*nCols+k+0]),c03v,sumi1v);
            sumi2v = fmadd256d(load256d012(&ddNTdtdNT[(i+2)*nCols+k+0]),c03v,sumi2v);
        }

        __m256d sumi03=  hsum4x4(sumi0v,sumi1v,sumi2v,sumi2v);
        __m256d dNi03dt = load256d012(&dNTdt[i+0]);
        __m256d result = sub256d(dNi03dt,sumi03);
        result = mul256d(VByNtotv,result);
        __m256d result0v = _mm256_permute4x64_pd(result,0b00000000);
        __m256d result1v = _mm256_permute4x64_pd(result,0b01010101);
        __m256d result2v = _mm256_permute4x64_pd(result,0b10101010);
        __m256d hai03v = load256d012(&ha[i]);

        __m256d hai0v = _mm256_set1_pd(ha[i+0]);
        __m256d hai1v = _mm256_set1_pd(ha[i+1]);
        __m256d hai2v = _mm256_set1_pd(ha[i+2]);
        for (int j=0; j<Ns-3; j=j+4)
        {
            __m256d Jaci0j03v = add256d(result0v,load256d(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d(&ddNTdtdNT[(i+2)*nCols+j+0]));

            store256d(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);
            store256d(&ddNTdtdNT[(i+2)*nCols+j+0],-Jaci2j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            store256d(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }
        {
            int j=Ns-3;
            __m256d Jaci0j03v = add256d(result0v,load256d012(&ddNTdtdNT[(i+0)*nCols+j+0]));
            __m256d Jaci1j03v = add256d(result1v,load256d012(&ddNTdtdNT[(i+1)*nCols+j+0]));
            __m256d Jaci2j03v = add256d(result2v,load256d012(&ddNTdtdNT[(i+2)*nCols+j+0]));

            store256d012(&ddNTdtdNT[(i+0)*nCols+j+0],-Jaci0j03v);
            store256d012(&ddNTdtdNT[(i+1)*nCols+j+0],-Jaci1j03v);
            store256d012(&ddNTdtdNT[(i+2)*nCols+j+0],-Jaci2j03v);

            __m256d ddTdtdNj =  fnmadd256d(hai0v,Jaci0j03v,load256d012(&ddNTdtdNT[Ns*nCols+j]));
            ddTdtdNj =          fnmadd256d(hai1v,Jaci1j03v,ddTdtdNj);
            ddTdtdNj =          fnmadd256d(hai2v,Jaci2j03v,ddTdtdNj);
            store256d012(&ddNTdtdNT[Ns*nCols+j],ddTdtdNj);
        }

        ddNTdtdNT[(i+0)*nCols+i+0] += scale;
        ddNTdtdNT[(i+1)*nCols+i+1] += scale;
        ddNTdtdNT[(i+2)*nCols+i+2] += scale;
        // Convert dcdt to dNidt
        __m256d dNi03dtv = load256d012(&dNTdt[i+0]);
        dNi03dtv = _mm256_mul_pd(Vv,dNi03dtv);
        store256d012(&dNTdt[i],dNi03dtv);

        __m256d ddNidtdTv = _mm256_setr_pd(ddNTdtdNT[(i+0)*nCols+Ns],ddNTdtdNT[(i+1)*nCols+Ns],ddNTdtdNT[(i+2)*nCols+Ns],0);
        ddNidtdTv = fnmadd256d(sumi03,_mm256_set1_pd(invT),ddNidtdTv);
        ddNidtdTv = mul256d(ddNidtdTv,_mm256_set1_pd(V));
        ddNidtdTv = fmadd256d(dNi03dtv,_mm256_set1_pd(invT),ddNidtdTv);
        ddNTdtdNT[(i+0)*nCols+Ns] = -get0(ddNidtdTv);
        ddNTdtdNT[(i+1)*nCols+Ns] = -get1(ddNidtdTv);
        ddNTdtdNT[(i+2)*nCols+Ns] = -get2(ddNidtdTv);

        ddTdtdTv = fnmadd256d(ddNidtdTv,hai03v,ddTdtdTv);
        ddTdtdTv = fnmadd256d(dNi03dtv,load256d012(&cp[i+0]),ddTdtdTv);
    }

    const double dTdt = dNTdt[Ns];
    dTdtv = _mm256_set1_pd(dNTdt[Ns]);
    __m256d invcptotv = _mm256_set1_pd(invcptot);
    for (int i=0; i<Ns-3; i=i+4)
    {
        __m256d ddTdtdNi0v = fnmadd256d(load256d(&cp[i+0]),dTdtv,load256d(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d(&ddNTdtdNT[Ns*nCols+i],-ddTdtdNi0v);
    }
    {
        int i=Ns-3;
        __m256d ddTdtdNi0v = fnmadd256d(load256d012(&cp[i+0]),dTdtv,load256d012(&ddNTdtdNT[Ns*nCols+i+0]));
        ddTdtdNi0v = mul256d(ddTdtdNi0v,invcptotv);
        store256d012(&ddNTdtdNT[Ns*nCols+i],-ddTdtdNi0v);
    }
    double ddTdtdT = hsum4(ddTdtdTv);
    ddTdtdT = ddTdtdT - dTdt*dcptotdT;
    ddTdtdT = ddTdtdT*invcptot;
    ddNTdtdNT[Ns*nCols+Ns] = -ddTdtdT + scale;
}

}// End of namespace Foam
