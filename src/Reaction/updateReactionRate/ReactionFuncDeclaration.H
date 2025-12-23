// Note: This file is for text substitution purposes only, 
// not compliant with standard C++ header file usage conventions

// These function is used for computing molar reaction rate
// RR:  reversible reaction
// IR:  irreversible reaction reverse rate constant is zero
// NER: non-equilibrium reaction, reverse rate constant is computed using
//      Arrhenius form instead of equilibrium rate constant
// GNI: non-integer reaction. e.g. 1.5A^1.5 = B
// GI:  integer reaction with number of reactants > 3 or 
//      number of products > 3. e.g. 3A+B = C
            
            void RFGNI
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RFGI
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF11RR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF11IR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF11NER
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF12RR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF12IR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF12NER
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF13RR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF13IR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF13NER
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF21RR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF21IR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF21NER
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF22RR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF22IR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF22NER
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF23RR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF23IR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF23NER
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF31RR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF31IR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF31NER
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF32RR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF32IR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF32NER
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void RF33RR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF33IR
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void RF33NER
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;