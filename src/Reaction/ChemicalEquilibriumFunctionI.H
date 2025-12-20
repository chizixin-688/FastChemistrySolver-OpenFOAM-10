inline void OptReaction::update_Pow_pByRT_SumVki
(
    const double T_
) const
{
    for(auto it = Pow_pByRT_SumVki_I.begin(); it!=Pow_pByRT_SumVki_I.end();it++)
    {
        int sumVki = it->first;
        double val = 0;
        double PByRT = this->Pstd*this->invRu*this->invT;
        {
            if (sumVki < 0)
            {
                PByRT = 1.0 / PByRT;
                sumVki = -sumVki;
            }

            double r = 1.0;
            while (sumVki)
            {
                if (sumVki & 1) r *= PByRT;
                PByRT *= PByRT;
                sumVki >>= 1;
            }
            val = r;
        }
        it->second = val;
    }
}

inline void OptReaction::update_Pow_pByRT_SumVki2
(
    const double T_
) const
{
    //const double PByRT = this->Pstd/(this->Ru*T_);
    const double PByRT = this->Pstd*this->invRu*this->invT;
    Pow_pByRT_SumVki[0] = 1/(PByRT*PByRT);
    Pow_pByRT_SumVki[1] = 1/(PByRT);      
    Pow_pByRT_SumVki[2] = 1;   
    Pow_pByRT_SumVki[3] = PByRT;
    Pow_pByRT_SumVki[4] = PByRT*PByRT;
}
