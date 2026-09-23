/*---------------------------------------------------------------------------*\
  Description
      Declaration of the molar-concentration-based reaction-rate update
      functions (RF = Reaction Function) of the gas-phase reaction object.
      This file is included in OptReaction.H for text substitution purposes
      only.

      Each function accumulates, for one class of elementary reactions, the
      contribution to the net rate of production dNdtByV
      (mol/(m^3*s)) using the precomputed forward rate constants Kf_.

      Naming convention:  RF <Nr> <Np> <Type>
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
          RFGI  : general reaction with integer stoichiometric coefficients
                  (more than 3 reactants and/or more than 3 products)
          RFGNI : general reaction with non-integer stoichiometric
                  coefficients, e.g. 1.5 A^1.5 = B

      Common parameters of every function:
        c              : molar concentrations of the species
        dNdtByV        : net rate of production (accumulated in-place)
        ExpNegGbyRT    : exp(-G/(RT)) of every species, used to build the
                         chemical equilibrium constant Kc

      The functions are dispatched through the RFptr member (populated in
      createGasPhaseReactionFromFoamDict.C) and invoked from
      OptReaction::ReactionRate() and OptReaction::Tc().

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

            // ----------------------------------------------------------------------
            // 1-1 reactions: one reactant, one product. e.g. A <=> B.
            // ----------------------------------------------------------------------

            // RF11RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void RF11RR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF11IR: Irreversible reaction; the reverse rate constant is zero.

            void RF11IR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF11NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void RF11NER
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // 1-2 reactions: one reactant, two products. e.g. A <=> B+C.
            // ----------------------------------------------------------------------

            // RF12RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void RF12RR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF12IR: Irreversible reaction; the reverse rate constant is zero.

            void RF12IR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF12NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void RF12NER
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // 1-3 reactions: one reactant, three products. e.g. A <=> B+C+D.
            // ----------------------------------------------------------------------

            // RF13RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void RF13RR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF13IR: Irreversible reaction; the reverse rate constant is zero.

            void RF13IR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF13NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void RF13NER
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // 2-1 reactions: two reactants, one product. e.g. A+B <=> C.
            // ----------------------------------------------------------------------

            // RF21RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void RF21RR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF21IR: Irreversible reaction; the reverse rate constant is zero.

            void RF21IR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF21NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void RF21NER
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // 2-2 reactions: two reactants, two products. e.g. A+B <=> C+D.
            // ----------------------------------------------------------------------

            // RF22RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void RF22RR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF22IR: Irreversible reaction; the reverse rate constant is zero.

            void RF22IR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF22NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void RF22NER
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // 2-3 reactions: two reactants, three products. e.g. A+B <=> C+D+E.
            // ----------------------------------------------------------------------

            // RF23RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void RF23RR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF23IR: Irreversible reaction; the reverse rate constant is zero.

            void RF23IR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF23NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void RF23NER
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // 3-1 reactions: three reactants, one product. e.g. A+B+C <=> D.
            // ----------------------------------------------------------------------

            // RF31RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void RF31RR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF31IR: Irreversible reaction; the reverse rate constant is zero.

            void RF31IR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF31NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void RF31NER
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // 3-2 reactions: three reactants, two products. e.g. A+B+C <=> D+E.
            // ----------------------------------------------------------------------

            // RF32RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void RF32RR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF32IR: Irreversible reaction; the reverse rate constant is zero.

            void RF32IR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF32NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void RF32NER
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // 3-3 reactions: three reactants, three products. e.g. A+B+C <=> D+E+F.
            // ----------------------------------------------------------------------

            // RF33RR: Reversible reaction; the reverse rate constant is obtained from the chemical equilibrium constant.

            void RF33RR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF33IR: Irreversible reaction; the reverse rate constant is zero.

            void RF33IR
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // RF33NER: Non-equilibrium reversible reaction; the reverse rate constant is computed with the Arrhenius form instead of the equilibrium constant.

            void RF33NER
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // ----------------------------------------------------------------------
            // General reactions with integer stoichiometric coefficients
            // (more than 3 reactants and/or more than 3 products).
            // e.g. A+A+A+A <=> B+C+D, A <=> D+E+F+G.
            // ----------------------------------------------------------------------

            // RFGI: general reaction with integer stoichiometric coefficients.
            // It handles the reversible / irreversible / non-equilibrium
            // reverse-rate cases internally according to isIrreversible[i].

            void RFGI
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

            // ----------------------------------------------------------------------
            // General reactions with non-integer stoichiometric coefficients.
            // e.g. 3.5 A^1.4 <=> B+C+D.
            // ----------------------------------------------------------------------

            // RFGNI: general reaction with non-integer stoichiometric coefficients.
            // It handles the reversible / irreversible / non-equilibrium
            // reverse-rate cases internally according to isIrreversible[i].

            void RFGNI
            (
                const double* __restrict__ c,
                double* __restrict__ dNdtByV,
                const double* __restrict__ ExpNegGbyRT
            )const noexcept;

