/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#include "FastChemistryModel.H"
#include "basicChemistryModel.H"
#include "OptReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::FastChemistryModel<ThermoType>::FastChemistryModel
(
    const fluidReactionThermo& thermo
)
:   basicChemistryModel(thermo),
    Yvf_(this->thermo().composition().Y()),
    nSpecie_(Yvf_.size()),
    jacobianType_
    (
        this->found("jacobian")
      ? jacobianTypeNames_.read(this->lookup("jacobian"))
      : jacobianType::exact
    ),
    reaction(),
    RR_(this->nSpecie_),
    n_(this->nSpecie_+1),
    //alignN(n_+(4-(this->nSpecie_+1)%4)),
    alignN(((n_+3)/4)*4),
    Treact(this->lookupOrDefault("Treact",0)),
    DLBthreshold(this->lookupOrDefault("DLBthreshold",1.0)),
    MaxIter(this->lookupOrDefault("Iter",1)),
    cpuLoadTransferTable(Pstream::nProcs()),
    CPUtimeField(thermo.T().mesh().C().size()),
    chemistryIntegrationTime(Pstream::nProcs()),
    sendBufferSize_(Pstream::nProcs()),
    recvBufferSize_(Pstream::nProcs()),
    sendBuffer_(Pstream::nProcs()),
    recvBuffer_(Pstream::nProcs()),
    recvBufPos_(Pstream::nProcs()),
    firstTime(true),
    skip(thermo.T().mesh().C().size(),false),
    IamBusyProcess(Pstream::nProcs(),true),
    Balance(this->lookupOrDefault("balance", false))
{

    const IOdictionary thermoDict
    (
        physicalProperties::findModelDict(this->mesh(), word::null)
    );
    const IOdictionary chemistryProperties
    (
        IOobject
        (
            "chemistryProperties",
            this->mesh().time().constant(),
            this->mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    label defaultIndex = 0;
    const word defaultSpecie = thermoDict.lookup("defaultSpecie");
    Info<<"The default specie is "<<defaultSpecie<<endl;

    hashedWordList speciesTable(thermoDict.lookup("species"));
    defaultIndex = speciesTable[defaultSpecie];

    if(defaultIndex<0 ||defaultIndex>=this->nSpecie())
    {
        FatalErrorInFunction
                    << "Index of default species is wrong!"
                    << Foam::abort(FatalError);
    }
    reaction.readInfo(chemistryProperties,thermoDict);
    FastChemistry::createIdealGasFromFoamDict(thermoDict,gas);
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField
            (
                IOobject
                (
                    "RR." + Yvf_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                thermo.T().mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "FastChemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction() << endl;

    {
        

        //size_t alignN = N;
        //size_t n_ = (this->nSpecie()+1);
        size_t totalSize = 12*alignN + 3*alignN*n_;
        size_t bytes = totalSize * sizeof(double);
        if (posix_memalign(reinterpret_cast<void**>(&this->buffer), 32, bytes))
        {
            throw std::bad_alloc();
        }
        std::memset(this->buffer, 0, bytes);
        size_t pos = 0;

        for (int i = 0; i < 12; i++)
        {
            YTpWork[i] = buffer + pos;
            pos   += alignN;
        }
        for (int i = 0; i < 3; i++)
        {
            YTpYTpWork[i] = buffer + pos;
            pos   += alignN * n_;
        }
    }
    forAll(cpuLoadTransferTable,i)
    {
        cpuLoadTransferTable[i].resize(Pstream::nProcs(),0);
    }

    reaction.alignN = this->alignN;

    // select the jacobian function
    {
        int remain = this->nSpecie()%8;

        if(remain==0)
        {
            if(jacobianType_==jacobianType::exact){this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::ddYdtdY_Vec844_0;}
            else{this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_0;}
            this->ddYdtdTPtr = &FastChemistryModel<ThermoType>::ddYdtdT_Vec88_0;
            this->ddTdtdYTPtr = &FastChemistryModel<ThermoType>::ddTdtdYT_Vec88_0;
        }
        else if(remain==1)
        {
            if(jacobianType_==jacobianType::exact){this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::ddYdtdY_Vec844_1;}
            else{this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_1;}
            this->ddYdtdTPtr = &FastChemistryModel<ThermoType>::ddYdtdT_Vec88_1;
            this->ddTdtdYTPtr = &FastChemistryModel<ThermoType>::ddTdtdYT_Vec88_1;
        }
        else if(remain==2)
        {
            if(jacobianType_==jacobianType::exact){this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::ddYdtdY_Vec844_2;}
            else{this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_2;}
            this->ddYdtdTPtr = &FastChemistryModel<ThermoType>::ddYdtdT_Vec88_2;
            this->ddTdtdYTPtr = &FastChemistryModel<ThermoType>::ddTdtdYT_Vec88_2;
        }
        else if(remain==3)
        {
            if(jacobianType_==jacobianType::exact){this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::ddYdtdY_Vec844_3;}
            else{this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_3;}
            this->ddYdtdTPtr = &FastChemistryModel<ThermoType>::ddYdtdT_Vec88_3;
            this->ddTdtdYTPtr = &FastChemistryModel<ThermoType>::ddTdtdYT_Vec88_3;
        }
        else if(remain==4)
        {
            if(jacobianType_==jacobianType::exact){this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::ddYdtdY_Vec844_4;}
            else{this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_0;}
            this->ddYdtdTPtr = &FastChemistryModel<ThermoType>::ddYdtdT_Vec88_4;
            this->ddTdtdYTPtr = &FastChemistryModel<ThermoType>::ddTdtdYT_Vec88_4;
        }
        else if(remain==5)
        {
            if(jacobianType_==jacobianType::exact){this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::ddYdtdY_Vec844_5;}
            else{this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_1;}
            this->ddYdtdTPtr = &FastChemistryModel<ThermoType>::ddYdtdT_Vec88_5;
            this->ddTdtdYTPtr = &FastChemistryModel<ThermoType>::ddTdtdYT_Vec88_5;
        }
        else if(remain==6)
        {
            if(jacobianType_==jacobianType::exact){this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::ddYdtdY_Vec844_6;}
            else{this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_2;}
            this->ddYdtdTPtr = &FastChemistryModel<ThermoType>::ddYdtdT_Vec88_6;
            this->ddTdtdYTPtr = &FastChemistryModel<ThermoType>::ddTdtdYT_Vec88_6;
        }
        else
        {
            if(jacobianType_==jacobianType::exact){this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::ddYdtdY_Vec844_7;}
            else{this->ddYdtdYPtr = &FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_3;}
            this->ddYdtdTPtr = &FastChemistryModel<ThermoType>::ddYdtdT_Vec88_7;
            this->ddTdtdYTPtr = &FastChemistryModel<ThermoType>::ddTdtdYT_Vec88_7;
        }
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::FastChemistryModel<ThermoType>::~FastChemistryModel()
{
    free(this->buffer);

    for (int i = 0; i < 12; i++)
    {
        YTpWork[i] = nullptr;
    }
    for (int i = 0; i < 3; i++)
    {
        YTpYTpWork[i] = nullptr;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class ThermoType>
void Foam::FastChemistryModel<ThermoType>::derivatives
(
    const scalar t,
    const label li,
    const double p,
    double* __restrict__ Phi,    
    double* __restrict__ dPhidt,
    double* __restrict__ Cp,
    double* __restrict__ Ha
) const
{
    double* __restrict__ c = YTpYTpWork[0];

    // Constrain temperature to valid range (given by thermo.dat)
    const double Tlowmin = gas->TlowMin;
    const double Thighmax = gas->ThighMax;
    double T = Phi[this->nSpecie()];
    T = T<Tlowmin?Tlowmin:T;
    T = T>Thighmax?Thighmax:T;

    // Computing temperature and pressure, required by thermo and chemistry
    __m256d TP = _mm256_setr_pd(T,p,1,1);
    TP = vec256_logd(TP);
    const double invT = 1.0/T;
    const double logT = get_elem0(TP);
    const double logP = get_elem1(TP);
    gas->invT = invT;
    gas->logT = logT;
    reaction.invT = invT;
    reaction.logT = logT;
    reaction.logP = logP;

    // set to zero
    std::memset(dPhidt, 0, alignN*sizeof(double));
    
    

    /*for (int i=0; i<this->nSpecie(); i++)
    {
        Phi[i] = std::max(Phi[i], 0.0);
    }
    
    double rhoM = 0;
    double RuTByP = reaction.Ru*T/p;
    __m256d RuTByPv = _mm256_set1_pd(RuTByP);
    __m256d rhoMv = _mm256_setzero_pd();
    for (int i=0; i<nSpecie_-remain; i=i+4)
    {
        __m256d YTpv = load256d(&Phi[i+0]);
        __m256d invWv = load256d(&reaction.invW[i+0]);
        rhoMv = _mm256_fmadd_pd(_mm256_mul_pd(YTpv,invWv),RuTByPv,rhoMv);
    }
    for(int i = nSpecie_-remain; i<nSpecie_;i++)
    {
        rhoM += Phi[i]*reaction.invW[i]*RuTByP;            
    }
    rhoM += reaction.hsum4(rhoMv);
    double invrhoM = rhoM;
    rhoM = 1/rhoM;

    for (label i=0; i<nSpecie_; i ++)
    {
        c[i] = rhoM*reaction.invW[i]*Phi[i];
    }*/


    // Compute the thermodynamic parameters required by chemical reaction rate

    // i-th species variables (size:nSpecie):
    // c[i]             :Molar concentration                                [kmol/m^3]
    // Cp[i]            :Specific heat capacity                             [J/kg/K]
    // Ha[i]            :Absolute enthalpy                                  [J/kg]
    // tmp_Exp[i]       :Dimensionless standard Gibbs energy(-Gstd/Ru/T)    [dimLess] 

    // Mixture variable:
    // gas->rhoM        :Mixture density                                    [kg/m^3]
    // gas->vM          :Mixture specific volume                            [m^3/kg]
    // Cp[nSpecie()]    :Mixture specific heat capacity                     [J/kg/K]
    gas->DerivativeThermoYT(T,p,Phi,c,Cp,Ha,reaction.tmp_Exp,reaction.negGstdByRT);


    // Compute the molar reaction rate

    // i-th species variables (size:nSpecie):
    // dPhidt[i]        :Molar reaction rate                                [kmol/m^3/s]
    reaction.dNdtByV(p,T,c,dPhidt);


    // Compute the change rate of mass fraction and temperature

    // i-th species variables (size:nSpecie):
    // dPhidt[i]        :Mass fraction change rate                          [1/s]

    // Mixture variable:
    // dTdt             :Temperature change rate                            [K/s]
    double dTdt = 0;
    double vm = gas->vM;
    const double* __restrict__ W = gas->W;
    __m256d dTdtv = _mm256_setzero_pd();
    __m256d invrhoMv = _mm256_set1_pd(vm);
    int remain = this->nSpecie()%4;
    for (label i=0; i<this->nSpecie()-remain; i=i+4)
    {
        __m256d Wv = load256d(&W[i]);
        __m256d dYTdtv = load256d(&dPhidt[i]);

        dYTdtv = _mm256_mul_pd(_mm256_mul_pd(Wv,invrhoMv),dYTdtv);
        store256d(&dPhidt[i],dYTdtv);

        __m256d Hav = load256d(&Ha[i]);
        dTdtv = _mm256_fmadd_pd(Hav,dYTdtv,dTdtv);
    }
    for(label i = nSpecie_-remain;i<nSpecie_;i++)
    {
        dPhidt[i] =dPhidt[i]*W[i]*vm;
        dTdt -= dPhidt[i]*Ha[i];
    }
    dTdt = dTdt -(reaction.hsum4(dTdtv));
    dTdt /= Cp[this->nSpecie()];
    dPhidt[this->nSpecie()] = dTdt;


}

template<class ThermoType>
void Foam::FastChemistryModel<ThermoType>::jacobian
(
    const scalar t,
    const label li,
    const double p,
    double* __restrict__ Phi,
    double* __restrict__ dPhidt,
    double* __restrict__ Jac
) const 
{

    for(int i = 0; i < this->nSpecie();i++)
    {
        Phi[i] = std::max(Phi[i], 0.0);
    }

    // Constrain temperature to valid range (given by thermo.dat)
    const double Tlowmin = gas->TlowMin;
    const double Thighmax = gas->ThighMax;
    double T = Phi[this->nSpecie()];
    T = T<Tlowmin?Tlowmin:T;
    T = T>Thighmax?Thighmax:T;

    // Computing temperature and pressure, required by thermo and chemistry
    __m256d TP = _mm256_setr_pd(T,p,1,1);
    TP = vec256_logd(TP);
    const double invT = 1.0/T;
    const double logT = get_elem0(TP);
    const double logP = get_elem1(TP);
    gas->invT = invT;
    gas->logT = logT;
    reaction.invT = invT;
    reaction.logT = logT;
    reaction.logP = logP;


    
    double* __restrict__ ddNdtByVdcT = YTpYTpWork[0];
    double* __restrict__ c           = YTpWork[3];
    double* __restrict__ dBdT        = YTpWork[4];
    double* __restrict__ dCpdT       = YTpWork[5];
    double* __restrict__ Cp          = YTpWork[6];
    double* __restrict__ Ha          = YTpWork[7];
    double* __restrict__ WiByrhoM    = YTpWork[8];
    double* __restrict__ rhoMByRhoi      = YTpWork[10];
    double* __restrict__ dcdY      = YTpYTpWork[2];
    {
        size_t size = alignN*(this->nSpecie()+1);
        std::memset(ddNdtByVdcT, 0, size * sizeof(double));
    }
    {
        size_t size = alignN;
        std::memset(dPhidt, 0, size * sizeof(double));
    }




    //for(int i =0;i<this->nSpecie()+1;i++)
    //{
    //    std::cout<<c[i]<<" "<<dBdT[i]<<" "<<dCpdT[i]<<" "<<Cp[i]<<" "<<Ha[i]<<" "<<rhoMByRhoi[i]<<" "<<WiByrhoM[i]<<" "<<std::exp(reaction.tmp_Exp[i])<<std::endl;
    //}
    //std::exit(0);
    
    /*reaction.ddNdtByVdcTp
    (
        p,
        T,
        Phi,
        c,
        dPhidt,
        dBdT,
        dCpdT,
        Cp,
        Ha,
        rhoMByRhoi,
        WiByrhoM,
        ddNdtByVdcT
    );*/


    // d(dPhidt)dPhi. Phi=[Y0 Y1 Y2 ... YNs T]



    // Compute the thermodynamic parameters required by chemical Jacobian matrix
    // 

    // i-th species variables (size:nSpecie):
    // c[i]             :Molar concentration                                [kmol/m^3]
    // Cp[i]            :Specific heat capacity                             [J/kg/K]
    // Ha[i]            :Absolute enthalpy                                  [J/kg]
    // tmp_Exp[i]       :Dimensionless standard Gibbs energy(-Gstd/Ru/T)    [dimLess] 
    // dBdT[i]          :Partial derivative of -Gstd/Ru w.r.t temperature   [dimLess]
    // dCpdT[i]         :Partial derivative of Cp w.r.t temperature         [J/kg/K/K]
    // rhoMByRhoi[i]    :Mixture density divided by species density         [dimLess]
    // WiByrhoM[i]      :Specific molar volume                              [m^3/kmol]

    // Mixture variable:
    // gas->rhoM        :Mixture density                                    [kg/m^3]
    // gas->vM          :Mixture specific volume                            [m^3/kg]
    // Cp[nSpecie()]    :Mixture specific heat capacity                     [J/kg/K]
    gas->JacobianThermoYT
    (
        p,
        T,
        Phi,
        c,
        reaction.tmp_Exp,
        reaction.negGstdByRT,
        dBdT,
        dCpdT,
        Cp,
        Ha,
        rhoMByRhoi,
        WiByrhoM
    );

    // Compute the molar based jacobian matrix d(dcTdt)dcT
    // cT = [c0, c1, ..., cNs, T]
    reaction.ddNdtByVdcTp
    (
        p,
        T,
        Phi,
        c,
        dPhidt,
        dBdT,
        ddNdtByVdcT
    );

    // Compute the mass fraction based jacobian matrix d(dYTdt)dYT
    // YT = [Y0, Y1, ..., YNs, T]
    double* __restrict__ invW = gas->invW;
    double rhoM = gas->rhoM;
    double alphav = gas->alphav();

    (this->*ddYdtdYPtr)(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
    (this->*ddYdtdTPtr)(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac,alphav);
    (this->*ddTdtdYTPtr)(Cp,dCpdT,Ha,dPhidt,Jac);

/*for(int i = 0; i < this->nSpecie()+1;i++)
{
    for(int j = 0; j < this->nSpecie()-1+1;j++)
    {
        std::cout<<Jac[i*(this->alignN)+j]<<" ";
    }
    int j = this->nSpecie()-1+1;
    std::cout<<Jac[i*(this->alignN)+j]<<std::endl;
}
std::exit(0);*/

    /*int remain = this->nSpecie()%4;
    switch (jacobianType_)
    {
        case jacobianType::fast:
            if(remain==0)
            {
                //reaction.FastddYdtdY_Vec0(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,Jac);
                this->FastddYdtdY_Vec44_0(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                reaction.ddYdtdTP_Vec_0(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac);
                reaction.ddTdtdYT_Vec_0(Cp,dCpdT,Ha,dPhidt,Jac);
            }
            else if(remain==1)
            {
                //reaction.FastddYdtdY_Vec1(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,Jac);
                this->FastddYdtdY_Vec44_1(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                reaction.ddYdtdTP_Vec_1(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac);  
                reaction.ddTdtdYT_Vec_1(Cp,dCpdT,Ha,dPhidt,Jac);
            }
            else if(remain==2)
            {
                //reaction.FastddYdtdY_Vec2(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,Jac);
                this->FastddYdtdY_Vec44_2(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                reaction.ddYdtdTP_Vec_2(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac);
                reaction.ddTdtdYT_Vec_2(Cp,dCpdT,Ha,dPhidt,Jac);
            }
            else
            {
                //reaction.FastddYdtdY_Vec3(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,Jac);
                this->FastddYdtdY_Vec44_3(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                reaction.ddYdtdTP_Vec_3(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac); 
                reaction.ddTdtdYT_Vec_3(Cp,dCpdT,Ha,dPhidt,Jac);
            }
        break;            
        case jacobianType::exact:
            if(this->nSpecie()%8==0)
            {
                this->ddYdtdY_Vec844_0(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                this->ddYdtdT_Vec88_0(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac,alphav);
                this->ddTdtdYT_Vec88_0(Cp,dCpdT,Ha,dPhidt,Jac);

                break;

            }
            else if(this->nSpecie()%8==1)
            {
                this->ddYdtdY_Vec844_1(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                this->ddYdtdT_Vec88_1(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac,alphav);
                this->ddTdtdYT_Vec88_1(Cp,dCpdT,Ha,dPhidt,Jac);
                break;
            }
            else if(this->nSpecie()%8==2)
            {
                this->ddYdtdY_Vec844_2(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                this->ddYdtdT_Vec88_2(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac,alphav);
                this->ddTdtdYT_Vec88_2(Cp,dCpdT,Ha,dPhidt,Jac);
                break;
            }
            else if(this->nSpecie()%8==3)
            {
                this->ddYdtdY_Vec844_3(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                this->ddYdtdT_Vec88_3(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac,alphav);
                this->ddTdtdYT_Vec88_3(Cp,dCpdT,Ha,dPhidt,Jac);
                break;
            }
            else if(this->nSpecie()%8==4)
            {
                this->ddYdtdY_Vec844_4(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                this->ddYdtdT_Vec88_4(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac,alphav);
                this->ddTdtdYT_Vec88_4(Cp,dCpdT,Ha,dPhidt,Jac);
                break;
            }
            else if(this->nSpecie()%8==5)
            {
                this->ddYdtdY_Vec844_5(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                this->ddYdtdT_Vec88_5(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac,alphav);
                this->ddTdtdYT_Vec88_5(Cp,dCpdT,Ha,dPhidt,Jac);
                break;
            }
            else if(this->nSpecie()%8==6)
            {
                this->ddYdtdY_Vec844_6(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                this->ddYdtdT_Vec88_6(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac,alphav);
                this->ddTdtdYT_Vec88_6(Cp,dCpdT,Ha,dPhidt,Jac);
                break;
            }
            else 
            {
                this->ddYdtdY_Vec844_7(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,invW,Phi,dcdY,Jac,rhoM);
                this->ddYdtdT_Vec88_7(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac,alphav);
                this->ddTdtdYT_Vec88_7(Cp,dCpdT,Ha,dPhidt,Jac);
                break;
            }
            if(remain==0)
            {
                reaction.ddYdtdY_Vec1_0(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,Phi,Jac);
                reaction.ddYdtdTP_Vec_0(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac);
                reaction.ddTdtdYT_Vec_0(Cp,dCpdT,Ha,dPhidt,Jac);
            }
            else if(remain==1)
            {
                reaction.ddYdtdY_Vec1_1(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,Phi,Jac);
                reaction.ddYdtdTP_Vec_1(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac);  
                reaction.ddTdtdYT_Vec_1(Cp,dCpdT,Ha,dPhidt,Jac);
            }
            else if(remain==2)
            {
                reaction.ddYdtdY_Vec1_2(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,Phi,Jac);
                reaction.ddYdtdTP_Vec_2(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac);
                reaction.ddTdtdYT_Vec_2(Cp,dCpdT,Ha,dPhidt,Jac);
            }
            else
            {
                reaction.ddYdtdY_Vec1_3(ddNdtByVdcT,rhoMByRhoi,WiByrhoM,dPhidt,Phi,Jac);
                reaction.ddYdtdTP_Vec_3(ddNdtByVdcT,WiByrhoM,c,dPhidt,Jac); 
                reaction.ddTdtdYT_Vec_3(Cp,dCpdT,Ha,dPhidt,Jac);
            }
        break;
    }*/
}




template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::FastChemistryModel<ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        volScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    scalarField& tc = ttc.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    double* __restrict__ C = this->YTpWork[0];

    if (this->chemistry_)
    {

        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            scalar Ti = T[celli];
            const scalar pi = p[celli];

            for (label i=0; i<nSpecie_; i++)
            {
                C[i] = rhoi*Yvf_[i][celli]*gas->invW[i];
            }

            // Constrain temperature to valid range (given by thermo.dat)
            const double Tlowmin = gas->TlowMin;
            const double Thighmax = gas->ThighMax;
            Ti = Ti<Tlowmin?Tlowmin:Ti;
            Ti = Ti>Thighmax?Thighmax:Ti;

            // Computing temperature and pressure, required by thermo and chemistry
            __m256d TP = _mm256_setr_pd(Ti,pi,1,1);
            TP = vec256_logd(TP);
            const double invT = 1.0/Ti;
            const double logT = get_elem0(TP);
            const double logP = get_elem1(TP);
            gas->invT = invT;
            gas->logT = logT;
            reaction.invT = invT;
            reaction.logT = logT;
            reaction.logP = logP;

            // Computing dimensionless standard Gibbs energy
            gas->negGstdByRT(Ti,reaction.tmp_Exp,reaction.negGstdByRT);

            // A reaction's rate scale is calculated as it's molar
            // production rate divided by the total number of moles in the
            // system.
            //
            // The system rate scale is the average of the reactions' rate
            // scales weighted by the reactions' molar production rates. This
            // weighting ensures that dominant reactions provide the largest
            // contribution to the system rate scale.
            //
            // The system time scale is then the reciprocal of the system rate
            // scale.
            //
            // Contributions from forward and reverse reaction rates are
            // handled independently and identically so that reversible
            // reactions produce the same result as the equivalent pair of
            // irreversible reactions.
            scalar sumW = 0, sumWRateByCTot = 0;
            reaction.Tc(celli,pi,Ti,C,sumW,sumWRateByCTot);

            double sumc = 0;
            for(int i = 0; i < this->nSpecie();i++)
            {
                sumc += C[i];
            }
            tc[celli] =
                sumWRateByCTot == 0 ? vGreat : sumW/sumWRateByCTot*sumc;
        }
    }
    ttc.ref().correctBoundaryConditions();
    return ttc;

}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::FastChemistryModel<ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {

        scalarField& Qdot = tQdot.ref();

        forAll(Yvf_, i)
        {
            forAll(Qdot, celli)
            {
                //const scalar hi = specieThermos_[i].Hf();
                const double hi = gas->Hf[i];
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::FastChemistryModel<ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{

    FatalErrorInFunction
                    << "This function is not supported and should not be used"
                    << Foam::abort(FatalError);

    tmp<volScalarField::Internal> tRR
    (
        volScalarField::Internal::New
        (
            "RR",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );
    return tRR;
}


template<class ThermoType>
void Foam::FastChemistryModel<ThermoType>::calculate()
{

    if (!this->chemistry_)
    {
        return;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();
    double* __restrict__ C = YTpWork[0];
    double* __restrict__ dNdtByV = YTpWork[1];
    //double* __restrict__ Cp = YTpWork[2];
    //double* __restrict__ Ha = YTpWork[3];

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        scalar Ti = T[celli];
        const scalar pi = p[celli];

        // Constrain temperature to valid range (given by thermo.dat)
        const double Tlowmin = gas->TlowMin;
        const double Thighmax = gas->ThighMax;
        Ti = Ti<Tlowmin?Tlowmin:Ti;
        Ti = Ti>Thighmax?Thighmax:Ti;

        // Computing temperature and pressure, required by thermo and chemistry
        __m256d TP = _mm256_setr_pd(Ti,pi,1,1);
        TP = vec256_logd(TP);
        const double invT = 1.0/Ti;
        const double logT = get_elem0(TP);
        const double logP = get_elem1(TP);
        gas->invT = invT;
        gas->logT = logT;
        reaction.invT = invT;
        reaction.logT = logT;
        reaction.logP = logP;

        // Computing molar concentration
        for (int i=0; i<this->nSpecie(); i++)
        {
            const scalar Yi = Yvf_[i][celli];
            C[i] = rhoi*Yi*gas->invW[i];
        }
        gas->negGstdByRT(Ti,reaction.tmp_Exp,reaction.negGstdByRT);

        std::memset(C, 0, this->alignN * sizeof(double)*4);

        reaction.dNdtByV(pi,Ti,C,dNdtByV);
        //reaction.dNdtByV(pi,Ti,C,dNdtByV,Cp,Ha);
        for (int i=0; i<this->nSpecie(); i++)
        {
            RR_[i][celli] = dNdtByV[i]*gas->W[i];
        }
    }
    return;
}


#include "FastChemistryModel_transientSolve.H"
#include "FastChemistryModel_localEulerSolve.H"
template<class ThermoType>
void Foam::FastChemistryModel<ThermoType>::exchange
(
    const UList<DynamicList<char>>& sendBufs,
    const List<std::streamsize>& recvSizes,
    List<DynamicList<char>>& recvBufs,
    const int tag,
    const label comm,
    const bool block
)
{
    if (!contiguous<char>())
    {
        FatalErrorInFunction
            << "Continuous data only." << sizeof(char) << Foam::abort(FatalError);
    }

    if (sendBufs.size() != UPstream::nProcs(comm))
    {
        FatalErrorInFunction
            << "Size of list " << sendBufs.size()
            << " does not equal the number of processors "
            << UPstream::nProcs(comm)
            << Foam::abort(FatalError);
    }

    recvBufs.setSize(sendBufs.size());

    if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        label startOfRequests = Pstream::nRequests();

        forAll(recvSizes, proci)
        {
            std::streamsize nRecv = recvSizes[proci]; 


            if (proci != Pstream::myProcNo(comm) && nRecv > 0)
            {

                recvBufs[proci].setSize(static_cast<Foam::label>(nRecv)); 
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    reinterpret_cast<char*>(recvBufs[proci].begin()),
                    nRecv*sizeof(char),
                    tag,
                    comm
                );
            }
        }

        forAll(sendBufs, proci)
        {
            if (proci != Pstream::myProcNo(comm) && sendBufs[proci].size() > 0)
            {

                if
                (
                   !UOPstream::write
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        reinterpret_cast<const char*>(sendBufs[proci].begin()),
                        sendBufs[proci].size()*sizeof(char),
                        tag,
                        comm
                    )
                )
                {
                    FatalErrorInFunction
                        << "Cannot send outgoing message. "
                        << "to:" << proci << " nBytes:"
                        << label(sendBufs[proci].size()*sizeof(char))
                        << Foam::abort(FatalError);
                }
            }
        }

        if (block)
        {
            Pstream::waitRequests(startOfRequests); 
        }
    }

    recvBufs[Pstream::myProcNo(comm)] = sendBufs[Pstream::myProcNo(comm)];
}


// ************************************************************************* //
