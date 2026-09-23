/*---------------------------------------------------------------------------*\
  Description
      Declarations of the Jacobian kernel member functions of
      FastChemistryModel.

      This header is included inside the class definition (see
      FastChemistryModel.H) to keep the (large) list of kernel signatures
      readable.  Each declaration must match the definition provided by the
      corresponding file under src/ChemistryModel/Jacobian/.

      The kernels are grouped by state basis:

        - mass-fraction (Y) basis:
            species block        FastddYdtdY_Vec44_r (fast, Ns % 4)
                                 ddYdtdY_Vec844_r  (exact, Ns % 8)
            temperature column   ddYdtdT_Vec88_r   (d(dY/dt)/dT)
            temperature row      ddTdtdYT_Vec88_r  (d(dT/dt)/d(Y,T))

        - molar-number (N) basis:
            scalar reference     ddNdtdNT, ddTdtdNT, ddNTdtdNT
            SIMD fused kernels   ddNTdtdNT0..3 (and ...WithMergeI,
                                 which also builds A = -J + scale*I)

      The trailing digit r is the remainder of the number of species
      (Ns % 4 for the 4-lane kernels, Ns % 8 for the 8-lane kernels) that
      selects the corresponding specialised implementation.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

// ==== Y basis: species block, FAST approximation (4 lanes, Ns % 4) ====

    void FastddYdtdY_Vec44_0
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void FastddYdtdY_Vec44_1
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void FastddYdtdY_Vec44_2
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void FastddYdtdY_Vec44_3
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

// ==== Y basis: species block, EXACT (8 lanes, Ns % 8) ====

    void ddYdtdY_Vec844_0
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void ddYdtdY_Vec844_1
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void ddYdtdY_Vec844_2
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void ddYdtdY_Vec844_3
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void ddYdtdY_Vec844_4
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void ddYdtdY_Vec844_5
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void ddYdtdY_Vec844_6
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

    void ddYdtdY_Vec844_7
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ rhoMByRhoi,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ dPhidt,
        const double* __restrict__ invW,
        const double* __restrict__ Phi,
        double* __restrict__ dcdY,
        double* __restrict__ Jac,
        double rhoM
    ) const noexcept;

// ==== Y basis: temperature column  d(dY_i/dt)/dT  (Ns % 8) ====

    void ddYdtdT_Vec88_0
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ c,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac,
        double alphavM
    ) const noexcept;

    void ddYdtdT_Vec88_1
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ c,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac,
        double alphavM
    ) const noexcept;

    void ddYdtdT_Vec88_2
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ c,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac,
        double alphavM
    ) const noexcept;

    void ddYdtdT_Vec88_3
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ c,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac,
        double alphavM
    ) const noexcept;

    void ddYdtdT_Vec88_4
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ c,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac,
        double alphavM
    ) const noexcept;

    void ddYdtdT_Vec88_5
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ c,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac,
        double alphavM
    ) const noexcept;

    void ddYdtdT_Vec88_6
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ c,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac,
        double alphavM
    ) const noexcept;

    void ddYdtdT_Vec88_7
    (
        const double* __restrict__ ddNdtByVdcT,
        const double* __restrict__ WiByrhoM,
        const double* __restrict__ c,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac,
        double alphavM
    ) const noexcept;

// ==== Y basis: temperature row  d(dT/dt)/d(Y,T)  (Ns % 8) ====

    void ddTdtdYT_Vec88_0
    (
        const double* __restrict__ Cp,
        const double* __restrict__ dCpdT,
        const double* __restrict__ Ha,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac
    ) const noexcept;

    void ddTdtdYT_Vec88_1
    (
        const double* __restrict__ Cp,
        const double* __restrict__ dCpdT,
        const double* __restrict__ Ha,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac
    ) const noexcept;

    void ddTdtdYT_Vec88_2
    (
        const double* __restrict__ Cp,
        const double* __restrict__ dCpdT,
        const double* __restrict__ Ha,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac
    ) const noexcept;

    void ddTdtdYT_Vec88_3
    (
        const double* __restrict__ Cp,
        const double* __restrict__ dCpdT,
        const double* __restrict__ Ha,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac
    ) const noexcept;

    void ddTdtdYT_Vec88_4
    (
        const double* __restrict__ Cp,
        const double* __restrict__ dCpdT,
        const double* __restrict__ Ha,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac
    ) const noexcept;

    void ddTdtdYT_Vec88_5
    (
        const double* __restrict__ Cp,
        const double* __restrict__ dCpdT,
        const double* __restrict__ Ha,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac
    ) const noexcept;

    void ddTdtdYT_Vec88_6
    (
        const double* __restrict__ Cp,
        const double* __restrict__ dCpdT,
        const double* __restrict__ Ha,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac
    ) const noexcept;

    void ddTdtdYT_Vec88_7
    (
        const double* __restrict__ Cp,
        const double* __restrict__ dCpdT,
        const double* __restrict__ Ha,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac
    ) const noexcept;

// ==== N basis: scalar reference kernels ====

    void ddNdtdNT
    (
        const double* __restrict__ Phi,
        const double* __restrict__ c,
        double* __restrict__ dPhidt,
        double* __restrict__ Jac,
        double invNtot,
        double invT,
        double V,
        int Ns,
        int nCols
    ) const noexcept;

    void ddTdtdNT
    (
        double* __restrict__ dPhidt,
        const double* __restrict__ Phi,
        const double* __restrict__ ha,
        const double* __restrict__ cp,
        double* __restrict__ Jac,
        double invcptot,
        double dcptotdT,
        int Ns,
        int cols
    ) const noexcept;

// ==== N basis: scalar fused kernel (species block + T column + T row) ====

    void ddNTdtdNT
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
    ) const noexcept;

// ==== N basis: SIMD fused kernels (Ns % 4); ...WithMergeI also merges the implicit shift A = -J + scale*I ====

    void ddNTdtdNT0
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
    ) const noexcept;

    void ddNTdtdNT0WithMergeI
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
    ) const noexcept;

    void ddNTdtdNT1
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
    ) const noexcept;

    void ddNTdtdNT1WithMergeI
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
    ) const noexcept;

    void ddNTdtdNT2
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
    ) const noexcept;

    void ddNTdtdNT3
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
    ) const noexcept;

    void ddNTdtdNT3WithMergeI
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
    ) const noexcept;

    void ddNTdtdNT2WithMergeI
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
    ) const noexcept;
