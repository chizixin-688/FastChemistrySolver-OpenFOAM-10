// Note: This file is for text substitution purposes only, 
// not compliant with standard C++ header file usage conventions

// These function is used for computing mass fraction based 
// jacobian matrix using molar concentration based matrix

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