// Note: This file is for text substitution purposes only, 
// not compliant with standard C++ header file usage conventions

// These function is used for computing molar concentration based matrix
// RR:  reversible reaction
// IR:  irreversible reaction reverse rate constant is zero
// NER: non-equilibrium reaction, reverse rate constant is computed using
//      Arrhenius form instead of equilibrium rate constant

    // for one-one reactions. e.g. A<=>B.
    void JF11RR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF11IR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF11NER
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    // for one-two reactions. e.g. A<=>B+B, A<=>B+C.
    void JF12RR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF12IR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF12NER
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    // for one-three reactions. e.g. A<=>B+B+C, A<=>B+C+D.
    void JF13RR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF13IR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF13NER
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    // for two-one reactions. e.g. A<=>B+B+C, A<=>B+C+D.
    void JF21RR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF21IR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF21NER
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    // for two-two reactions. e.g. A+D<=>B+B, A+A<=>B+C.

    void JF22RR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF22IR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF22NER
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    // for two-three reactions. e.g. A+D<=>B+B, A+A<=>B+C.

    void JF23RR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF23IR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF23NER
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    // for three-one reactions. e.g. A+A+A<=>B, A+B+C<=>D.
    void JF31RR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF31IR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF31NER
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    // for three-two reactions. e.g. A+A+A<=>B+C, A+B+C<=>D+E.
    void JF32RR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF32IR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF32NER
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    // for three-three reactions. e.g. A+A+A<=>B+C+D, A+B+C<=>D+E+F.

    void JF33RR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF33IR
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    void JF33NER
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;


    // for general reactions with integer number. e.g. A+A+A+A<=>B+C+D, A<=>D+E+F+G.

    void updateJacobianGlobalIntegerReaction
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;

    // for general reactions with non-integer number. e.g. 3.5A^1.4<=>B+C+D.
    void updateJacobianGlobalNonIntegerReaction
    (
        const double* __restrict__ C,
        double* __restrict__ dNdtByV,
        double* __restrict__ ddNdtByVdcTp,
        const double* __restrict__ ExpNegGbyRT,
        const double* __restrict__ dBdT
    )const noexcept;