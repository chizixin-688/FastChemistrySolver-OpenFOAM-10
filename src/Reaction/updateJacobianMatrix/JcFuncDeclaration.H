/*---------------------------------------------------------------------------*\
  Description
              Declaration of the molar-concentration-based Jacobian update functions
              (JF = Jacobian Function) of the gas-phase reaction object. This file
              is included in OptReaction.H for text substitution purposes only.

              Each function accumulates, for one class of elementary reactions, the
              contributions to
                - the net rate of production dNdtByV, and
                - the molar-concentration-based Jacobian matrix ddNdtByVdcTp
                  (one row per species, alignN columns: nSpecies concentration
                  columns plus one extra column for the temperature derivative).

              Naming convention:  JF <Nr> <Np> <Type>
                Nr   : number of reactants  (1, 2 or 3)
                Np   : number of products   (1, 2 or 3)
                Type :
                  RR  : reversible reaction; the reverse rate constant is
                obtained from the chemical equilibrium constant
                  IR  : irreversible reaction; the reverse rate constant is zero
                  NER : non-equilibrium reaction; the reverse rate constant is
                computed with the Arrhenius form instead of the
                equilibrium constant
                In addition:
                  JFGI  : general reaction with integer stoichiometric coefficients
                  (more than 3 reactants and/or more than 3 products)
                  JFGNI : general reaction with non-integer stoichiometric
                  coefficients, e.g. 1.5 A^1.5 = B

              Common parameters of every function:
                C              : molar concentrations of the species
                dNdtByV        : net rate of production (accumulated in-place)
                ddNdtByVdcTp   : Jacobian matrix (accumulated in-place)
                ExpNegGbyRT    : exp(-G/(RuT)) of every species, used to build the
                         chemical equilibrium constant Kc
                dBdT           : temperature derivative of the standard-state Gibbs
                         energy terms, used for the temperature column of
                         the Jacobian

              The functions are dispatched through the JFptr member (populated in
              createGasPhaseReactionFromFoamDict.C) and invoked from
              OptReaction::ddNdtByVdcTp().

  Author
              Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

            // ----------------------------------------------------------------------
            // 1-1 reactions: one reactant, one product. e.g. A <=> B.
            // ----------------------------------------------------------------------

            // JF11RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void JF11RR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF11IR: Irreversible reaction; the reverse rate constant is zero.

            void JF11IR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF11NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void JF11NER
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // 1-2 reactions: one reactant, two products. e.g. A <=> B+C.
            // ----------------------------------------------------------------------

            // JF12RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void JF12RR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF12IR: Irreversible reaction; the reverse rate constant is zero.

            void JF12IR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF12NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void JF12NER
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // 1-3 reactions: one reactant, three products. e.g. A <=> B+C+D.
            // ----------------------------------------------------------------------

            // JF13RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void JF13RR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF13IR: Irreversible reaction; the reverse rate constant is zero.

            void JF13IR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF13NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void JF13NER
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // 2-1 reactions: two reactants, one product. e.g. A+B <=> C.
            // ----------------------------------------------------------------------

            // JF21RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void JF21RR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF21IR: Irreversible reaction; the reverse rate constant is zero.

            void JF21IR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF21NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void JF21NER
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // 2-2 reactions: two reactants, two products. e.g. A+B <=> C+D.
            // ----------------------------------------------------------------------

            // JF22RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void JF22RR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF22IR: Irreversible reaction; the reverse rate constant is zero.

            void JF22IR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF22NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void JF22NER
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // 2-3 reactions: two reactants, three products. e.g. A+B <=> C+D+E.
            // ----------------------------------------------------------------------

            // JF23RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void JF23RR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF23IR: Irreversible reaction; the reverse rate constant is zero.

            void JF23IR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF23NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void JF23NER
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // 3-1 reactions: three reactants, one product. e.g. A+B+C <=> D.
            // ----------------------------------------------------------------------

            // JF31RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void JF31RR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF31IR: Irreversible reaction; the reverse rate constant is zero.

            void JF31IR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF31NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void JF31NER
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // 3-2 reactions: three reactants, two products. e.g. A+B+C <=> D+E.
            // ----------------------------------------------------------------------

            // JF32RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void JF32RR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF32IR: Irreversible reaction; the reverse rate constant is zero.

            void JF32IR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF32NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void JF32NER
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // 3-3 reactions: three reactants, three products. e.g. A+B+C <=> D+E+F.
            // ----------------------------------------------------------------------

            // JF33RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void JF33RR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF33IR: Irreversible reaction; the reverse rate constant is zero.

            void JF33IR
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // JF33NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void JF33NER
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // ----------------------------------------------------------------------
            // General reactions with integer stoichiometric coefficients
            // (more than 3 reactants and/or more than 3 products).
            // e.g. A+A+A+A <=> B+C+D, A <=> D+E+F+G.
            // ----------------------------------------------------------------------

            // JFGI: general reaction with integer stoichiometric coefficients
            // (more than 3 reactants and/or more than 3 products). It handles the
            // reversible / irreversible / non-equilibrium reverse-rate cases
            // internally according to isIrreversible[i].

            void JFGI
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            // ----------------------------------------------------------------------
            // General reactions with non-integer stoichiometric coefficients.
            // e.g. 3.5 A^1.4 <=> B+C+D.
            // ----------------------------------------------------------------------

            // JFGNI: general reaction with non-integer stoichiometric coefficients.
            // It handles the reversible / irreversible / non-equilibrium reverse-rate
            // cases internally according to isIrreversible[i].

            void JFGNI
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

