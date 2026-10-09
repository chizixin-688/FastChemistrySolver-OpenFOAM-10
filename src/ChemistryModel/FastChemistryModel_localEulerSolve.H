/*---------------------------------------------------------------------------*\
  Description
      This file is for text substitution purposes only.

      Solving chemistry using local Euler time step

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

template<class UnusedThermo>
Foam::scalar Foam::FastChemistryModel<UnusedThermo>::solve
(
    const scalarField& deltaT
)
{

    auto GlobalTimeStart = std::chrono::high_resolution_clock::now();

    basicChemistryModel::correct();

    scalar deltaTMin = great;
    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho0vf(this->thermo().rho0());
    const volScalarField& rho0vf = trho0vf();

    const volScalarField& T0vf = this->thermo().T().oldTime();
    const volScalarField& p0vf = this->thermo().p().oldTime();

    double* __restrict__ Phi00 = this->YTpWork[0];
    double* __restrict__ Phi0  = this->YTpWork[1];

    if(firstTime || !Pstream::parRun() || !Balance)
    {
        firstTime = false;

        int64_t time = 0;

        forAll(rho0vf, celli)
        {
            auto TimeStart = std::chrono::high_resolution_clock::now();
            scalar T = T0vf[celli];
            if(T>this->Treact )
            {
                scalar p = p0vf[celli];
                const scalar rho0 = rho0vf[celli];

                // Load thermophysical variable into Phi vector
                for (label i=0; i<nSpecie_; i++)
                {
                    Phi00[i] = Yvf_[i].oldTime()[celli];
                }
                for (label i=0; i<nSpecie_; i++)
                {
                    Phi00[i] = max(0,Phi00[i]);
                    Phi0[i] =  Phi00[i];
                }
                Phi00[this->nSpecie()] = T;
                Phi0[this->nSpecie()]  = T;

                // Initialise time progress
                scalar timeLeft = deltaT[celli];

                // Calculate the chemical source terms
                while (timeLeft > small)
                {
                    scalar dt = timeLeft;
                    solve(celli, p, dt, deltaTChem_[celli]);
                    timeLeft -= dt;
                }

                deltaTMin = min(deltaTChem_[celli], deltaTMin);
                deltaTChem_[celli] = min(deltaTChem_[celli], deltaTChemMax_);


                // Set the RR vector (used in the solver)
                for (label i=0; i<nSpecie_; i++)
                {
                    Phi0[i] = max(0,Phi0[i]);
                }                
                for (label i=0; i<nSpecie_; i++)
                {
                    RR_[i][celli] = (Phi0[i]- Phi00[i])*rho0/deltaT[celli];
                }
            }
            else
            {
                for (label i=0; i<nSpecie_; i++)
                {
                    RR_[i][celli] = 0;
                }         
            }

            auto duration = (std::chrono::duration_cast<std::chrono::microseconds>
                (std::chrono::high_resolution_clock::now()-TimeStart));
            CPUtimeField[celli].first = static_cast<int64_t>(duration.count());
            CPUtimeField[celli].second = celli;
            time = time + static_cast<int64_t>(duration.count());
        }
        if(Pstream::parRun())
        {
            chemistryIntegrationTime[Pstream::myProcNo()] = time;
            Pstream::gatherList(chemistryIntegrationTime);
            if (Pstream::myProcNo()==0)
            {
                auto maxCPUtime = max(chemistryIntegrationTime);
                auto minCPUtime = min(chemistryIntegrationTime);    
                Info<<"Max/Min Unbalanced chemistry integration time: "<<maxCPUtime<<"/"<<minCPUtime<<endl;                    
            }
        }
        else
        {
            Info<<"Chemistry integration time: "<<time<<endl;
        }
        return deltaTMin;  
    }

    forAll(this->skip,i)
    {
        this->skip[i] = false;
    }

    List<label> busyProcs; 
    List<label> idleProcs;

    if(CurrentIter==0)
    {
        cpuLoadTransferTable.clear();        
        if (Pstream::myProcNo()==0)
        {

            int64_t averageCPUtime = 0;
            forAll(chemistryIntegrationTime,i)
            {
                averageCPUtime = averageCPUtime + chemistryIntegrationTime[i];
            }
            averageCPUtime = averageCPUtime/Pstream::nProcs() + 1;

            for(int i = 0; i < Pstream::nProcs(); i++)
            {
                if(
                    static_cast<double>(chemistryIntegrationTime[i])>
                    static_cast<double>(averageCPUtime)*DLBthreshold
                )
                {
                    busyProcs.append(i);
                    this->IamBusyProcess[i]=true;
                }
                else
                {
                    idleProcs.append(i);
                    this->IamBusyProcess[i]=false;
                }
            }

            List<int64_t> availableCPUtime(idleProcs.size());
            for(int i = 0; i < idleProcs.size(); i ++)
            {
                int j = idleProcs[i];
                availableCPUtime[i] = averageCPUtime - chemistryIntegrationTime[j];    
            }

            cpuLoadTransferTable.resize(Pstream::nProcs());
            forAll(cpuLoadTransferTable,i)
            {
                cpuLoadTransferTable[i].resize(Pstream::nProcs(),0);
            }

            for (int i = 0; i < busyProcs.size(); i++ )
            {
                int A = busyProcs[i];
                auto excessCPUtime = chemistryIntegrationTime[A] - averageCPUtime;
                for (int j = 0; j < idleProcs.size(); j++ )
                {
                    int B = idleProcs[j];

                    if(excessCPUtime == 0)
                    {
                        break;
                    }
                    else if(availableCPUtime[j]<=0)
                    {
                        continue;
                    }
                    auto CPUtimeToAssign = std::min(excessCPUtime,availableCPUtime[j]);

                    cpuLoadTransferTable[A][B] = CPUtimeToAssign;

                    excessCPUtime = excessCPUtime - CPUtimeToAssign;

                    availableCPUtime[j] = availableCPUtime[j] - CPUtimeToAssign;
                }
            }
        }
        Pstream::scatter(cpuLoadTransferTable);
        Pstream::scatter(IamBusyProcess);
    }

    CurrentIter++;
    if (CurrentIter >= MaxIter-1)
    {
        CurrentIter = 0;
    }

    List<int64_t> targetCPUtime(0);
    List<label> targetWorker(0);

    List<label> whoSendToMe(0);

    if(IamBusyProcess[Pstream::myProcNo()]==true)
    {
        const auto& myTable = cpuLoadTransferTable[Pstream::myProcNo()];
        forAll(myTable,i)
        {
            if (myTable[i]!=0)
            {
                targetWorker.append(i);
                targetCPUtime.append(myTable[i]);
            }
        }
    }
    if(IamBusyProcess[Pstream::myProcNo()]==false)
    {
        for(int i = 0; i < cpuLoadTransferTable.size(); i++)
        {
            if (cpuLoadTransferTable[i][Pstream::myProcNo()]!=0)
            {
                whoSendToMe.append(i);
            }
        }
    }

    int endPoint = 0;

    std::sort(CPUtimeField.begin(), CPUtimeField.end());
    std::reverse(CPUtimeField.begin(), CPUtimeField.end());
    

    List<List<dataBlock>> cellsToSend;        
    List<List<dataBlock>> cellsToReceive;     

    cellsToSend.clear();
    cellsToSend.resize( targetWorker.size() );

    if(IamBusyProcess[Pstream::myProcNo()]==true)
    {
        label startOfCell = 0;

        for(int i = 0; i < targetWorker.size(); i ++)
        {

            int64_t Timeleft = targetCPUtime[i];
            for(size_t j = startOfCell; j < CPUtimeField.size(); j ++) 
            {
                Timeleft = Timeleft - CPUtimeField[j].first;

                if (Timeleft > 0)
                {
                    int celli = CPUtimeField[j].second;
                    startOfCell++;
                    endPoint = endPoint + 1;
                    this->skip[celli]=true;
                    dataBlock data(nSpecie_);


                    data.T = T0vf[celli];
                    data.p = p0vf[celli];
                    for (label k=0; k<nSpecie_; k++)
                    {
                        data.Y_[k] =Yvf_[k].oldTime()[celli];
                    }

                    data.celli = celli;
                    data.deltaT = deltaT[celli];
                    data.deltaTChem_ = deltaTChem_[celli];

                    cellsToSend[i].append(data);
                }
                else
                {
                    break;
                }
            }
        }
    }
    cellsToReceive.clear();

    //********************************* MPI Communication ***********************************//
    {
        // Reset the buffer
        forAll(sendBufferSize_,i)
        {
            this->sendBufferSize_[i]=0;
            this->recvBufferSize_[i]=0;
            this->sendBuffer_[i].setSize(0);
            this->recvBuffer_[i].setSize(0);
            this->recvBufPos_[i] = 0;
        }

        label startOfRequests = Pstream::nRequests();

        if(IamBusyProcess[Pstream::myProcNo()]==true)
        {
            forAll(targetWorker,i)
            {

                // Create a stream for MPI communication
                UOPstream stream(
                                    Pstream::commsTypes::nonBlocking, 
                                    targetWorker[i], 
                                    sendBuffer_[targetWorker[i]],
                                    UPstream::msgType(),
                                    UPstream::worldComm,
                                    false//must be false to avoid send at destruct
                                );

                stream << cellsToSend[i];
                
                sendBufferSize_[i]=(sendBuffer_[targetWorker[i]].byteSize());

                if(
                    !UOPstream::write
                    (
                        UPstream::commsTypes::nonBlocking,
                        targetWorker[i],
                        reinterpret_cast<const char*>(&sendBufferSize_[i]),
                        sizeof(std::streamsize),
                        UPstream::msgType(),
                        UPstream::worldComm
                    )
                ){  FatalErrorInFunction
                    << "UOPstream::write failed!"
                    << Foam::abort(FatalError);}
            }
        }
        if(IamBusyProcess[Pstream::myProcNo()]==false)
        {
            forAll(whoSendToMe,i)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    whoSendToMe[i],
                    reinterpret_cast<char*>(&recvBufferSize_[whoSendToMe[i]]),
                    sizeof(std::streamsize),
                    UPstream::msgType(),
                    UPstream::worldComm
                );
            }
        }

        Pstream::waitRequests(startOfRequests);

        this->exchange
        (
            sendBuffer_,
            recvBufferSize_,
            recvBuffer_,
            UPstream::msgType(),//tag,
            UPstream::worldComm,//common,
            true //block=true
        );

        if(IamBusyProcess[Pstream::myProcNo()]==false)
        {
            forAll(whoSendToMe,i)
            {
                UIPstream stream(
                                    Pstream::commsTypes::nonBlocking,
                                    whoSendToMe[i], 
                                    recvBuffer_[whoSendToMe[i]],
                                    recvBufPos_[whoSendToMe[i]],
                                    UPstream::msgType(),
                                    UPstream::worldComm,
                                    true
                                );
                List<dataBlock> a_block(stream);
                cellsToReceive.append(a_block);
            }
        }
    }
    //********************************* MPI Communication ***********************************//

    if(IamBusyProcess[Pstream::myProcNo()]==false)
    {
        forAll(rho0vf, celli)
        {
            auto TimeStart = std::chrono::high_resolution_clock::now();
            scalar T = T0vf[celli];

            if(T>this->Treact)
            {
                const scalar rho0 = rho0vf[celli];
                scalar p = p0vf[celli];

                // Load thermophysical variable into Phi vector
                for (label i=0; i<nSpecie_; i++)
                {
                    Phi00[i] = Yvf_[i].oldTime()[celli];
                }
                for (label i=0; i<nSpecie_; i++)
                {
                    Phi00[i] = max(0,Phi00[i]);
                    Phi0[i] =  Phi00[i];
                }
                Phi00[this->nSpecie()] = T;
                Phi0[this->nSpecie()]  = T;

                // Initialise time progress
                scalar timeLeft = deltaT[celli];

                // Calculate the chemical source terms
                while (timeLeft > small)
                {
                    scalar dt = timeLeft;
                    {
                        solve(celli, p, dt, deltaTChem_[celli]);
                    }
                    timeLeft -= dt;
                }

                deltaTMin = min(deltaTChem_[celli], deltaTMin);
                deltaTChem_[celli] = min(deltaTChem_[celli], deltaTChemMax_);


                // Set the RR vector (used in the solver)
                for (label i=0; i<nSpecie_; i++)
                {
                    Phi0[i] = max(0,Phi0[i]);
                }
                for (label i=0; i<nSpecie_; i++)
                {
                    RR_[i][celli] = (Phi0[i]- Phi00[i])*rho0/deltaT[celli];
                }
            }
            else
            {
                for (label i=0; i<nSpecie_; i++)
                {
                    RR_[i][celli] = 0;
                }
            }

            auto duration = (std::chrono::duration_cast<std::chrono::microseconds>
                (std::chrono::high_resolution_clock::now()-TimeStart));
            CPUtimeField[celli].first = static_cast<int64_t>(duration.count());
            CPUtimeField[celli].second = celli;
        }
        forAll(cellsToReceive, A) 
        {
            forAll(cellsToReceive[A],celli)
            {
                auto TimeStart = std::chrono::high_resolution_clock::now();
                scalar T = cellsToReceive[A][celli].T;

                if(T>this->Treact)
                {
                    scalar p = cellsToReceive[A][celli].p;

                    // Load thermophysical variable into Phi vector
                    for (label i=0; i<nSpecie_; i++)
                    {
                        //Phi00[i] = Yvf_[i].oldTime()[celli];
                        Phi00[i] = cellsToReceive[A][celli].Y_[i];
                    }
                    for (label i=0; i<nSpecie_; i++)
                    {
                        Phi00[i] = max(0,Phi00[i]);
                        Phi0[i] =  Phi00[i];
                    }
                    Phi00[this->nSpecie()] = T;
                    Phi0[this->nSpecie()]  = T;

                    // Initialise time progress
                    scalar timeLeft = cellsToReceive[A][celli].deltaT;
                    scalar deltaTChem_celli = cellsToReceive[A][celli].deltaTChem_;

                    // Calculate the chemical source terms
                    while (timeLeft > small)
                    {
                        scalar dt = timeLeft;
                        
                        solve(celli, p, dt, deltaTChem_[celli]);
                        
                        timeLeft -= dt;
                    }
                    deltaTChem_celli = min(deltaTChem_celli, deltaTChemMax_);

                    for (label i=0; i<nSpecie_; i++)
                    {
                        cellsToReceive[A][celli].Y_[i] = max(0,Phi0[i]);
                    }
                }

                auto duration = (std::chrono::duration_cast<std::chrono::microseconds>
                    (std::chrono::high_resolution_clock::now()-TimeStart));

                cellsToReceive[A][celli].CPUtime = static_cast<int64_t>(duration.count());
            }
        }
    }
    if(IamBusyProcess[Pstream::myProcNo()]==true) 
    {
        forAll(rho0vf,celli)
        {
            if(this->skip[celli]==true)
            {
                continue;
            }

            auto TimeStart = std::chrono::high_resolution_clock::now();

            scalar T = T0vf[celli];

            if(T>this->Treact)
            {
                const scalar rho0 = rho0vf[celli];
                scalar p = p0vf[celli];

                // Load thermophysical variable into Phi vector
                for (label i=0; i<nSpecie_; i++)
                {
                    Phi00[i] = Yvf_[i].oldTime()[celli];
                }
                for (label i=0; i<nSpecie_; i++)
                {
                    Phi00[i] = max(0,Phi00[i]);
                    Phi0[i] =  Phi00[i];
                }
                Phi00[this->nSpecie()] = T;
                Phi0[this->nSpecie()]  = T;

                // Initialise time progress
                scalar timeLeft = deltaT[celli];

                // Calculate the chemical source terms
                while (timeLeft > small)
                {
                    scalar dt = timeLeft;
                    {
                        solve(celli, p, dt, deltaTChem_[celli]);
                    }
                    timeLeft -= dt;
                }

                deltaTMin = min(deltaTChem_[celli], deltaTMin);

                deltaTChem_[celli] = min(deltaTChem_[celli], deltaTChemMax_);

                // Set the RR vector (used in the solver)
                for (label i=0; i<nSpecie_; i++)
                {
                    Phi0[i] = max(0,Phi0[i]);
                }
                for (label i=0; i<nSpecie_; i++)
                {
                    RR_[i][celli] = (Phi0[i]- Phi00[i])*rho0/deltaT[celli];
                }
            }
            else
            {
                for (label i=0; i<nSpecie_; i++)
                {
                    RR_[i][celli] = 0;
                }
            }
            auto duration = (std::chrono::duration_cast<std::chrono::microseconds>
                (std::chrono::high_resolution_clock::now()-TimeStart));
            CPUtimeField[celli].first = static_cast<int64_t>(duration.count());
            CPUtimeField[celli].second = celli;
        }
    }

    //********************************* MPI Communication ***********************************//
    {
        forAll(sendBufferSize_,i)
        {
            this->sendBufferSize_[i]=0;
            this->recvBufferSize_[i]=0;
            this->sendBuffer_[i].setSize(0);
            this->recvBuffer_[i].setSize(0);
            this->recvBufPos_[i] = 0;
        }

        label startOfRequests = Pstream::nRequests();

        if(IamBusyProcess[Pstream::myProcNo()]==false)
        {
            forAll(whoSendToMe,i)
            {

                // Create a stream for MPI communication
                UOPstream stream(
                                    Pstream::commsTypes::nonBlocking, 
                                    whoSendToMe[i], 
                                    sendBuffer_[whoSendToMe[i]],
                                    UPstream::msgType(),
                                    UPstream::worldComm,
                                    false//must be false to avoid send at destruct
                                );
                stream << cellsToReceive[i];
                sendBufferSize_[i]=(sendBuffer_[whoSendToMe[i]].byteSize());

                if(
                    !UOPstream::write
                    (
                        UPstream::commsTypes::nonBlocking,
                        whoSendToMe[i],
                        reinterpret_cast<const char*>(&sendBufferSize_[i]),
                        sizeof(std::streamsize),
                        UPstream::msgType(),
                        UPstream::worldComm
                    )
                ){  FatalErrorInFunction
                    << "UOPstream::write failed!"
                    << Foam::abort(FatalError);}
            }
        }

        if(IamBusyProcess[Pstream::myProcNo()]==true) 
        {
            forAll(targetWorker,i)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    targetWorker[i],
                    reinterpret_cast<char*>(&recvBufferSize_[targetWorker[i]]),
                    sizeof(std::streamsize),
                    UPstream::msgType(),
                    UPstream::worldComm
                );
            }

        }
        Pstream::waitRequests(startOfRequests);

        // Send the data
        this->exchange
        (
            sendBuffer_,
            recvBufferSize_,
            recvBuffer_,
            UPstream::msgType(),//tag,
            UPstream::worldComm,//common,
            true //block=true
        );

        if(IamBusyProcess[Pstream::myProcNo()]==true)
        {

            cellsToSend.clear();
            forAll(targetWorker,i)
            {
                UIPstream stream(
                                Pstream::commsTypes::nonBlocking,
                                targetWorker[i], 
                                recvBuffer_[targetWorker[i]],
                                recvBufPos_[targetWorker[i]],
                                UPstream::msgType(),
                                UPstream::worldComm,
                                true
                                );
                List<dataBlock> a_block(stream);
                cellsToSend.append(a_block);
            }
        }
    }
    //********************************* MPI Communication ***********************************//


    if(IamBusyProcess[Pstream::myProcNo()]==true)
    {
        forAll(cellsToSend,A)
        {
            forAll(cellsToSend[A],j)
            {
                label celli = cellsToSend[A][j].celli;
                deltaTChem_[celli] = cellsToSend[A][j].deltaTChem_;

                deltaTMin = min(deltaTChem_[celli], deltaTMin);

                // Set the RR vector (used in the solver)
                for (label i=0; i<nSpecie_; i++)
                {
                    RR_[i][celli] = (cellsToSend[A][j].Y_[i]-Yvf_[i].oldTime()[celli])*rho0vf[celli]/deltaT[celli];                    
                }

                CPUtimeField[celli].first = cellsToSend[A][j].CPUtime;
                CPUtimeField[celli].second = celli;
            }
        }
    }

    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>
        (std::chrono::high_resolution_clock::now()-GlobalTimeStart));
    List<double> cpuTimesPerProcessor(Pstream::nProcs());

    cpuTimesPerProcessor[Pstream::myProcNo()] = static_cast<double>(duration.count());
    Pstream::gatherList(cpuTimesPerProcessor);

    int64_t totalCPUtime = 0;
    for(size_t i = 0; i < CPUtimeField.size(); i++ )
    {
        totalCPUtime = totalCPUtime + CPUtimeField[i].first;
    }
    chemistryIntegrationTime[Pstream::myProcNo()] = totalCPUtime;
    Pstream::gatherList(chemistryIntegrationTime);

    if(Pstream::myProcNo()==0)
    {
        auto maxUnbalancedCPUtime = max(chemistryIntegrationTime);
        auto minUnbalancedCPUtime = min(chemistryIntegrationTime); 
        
        auto maxbalancedCPUtime = max(cpuTimesPerProcessor);
        auto minbalancedCPUtime = min(cpuTimesPerProcessor);       
        Info<<"Max/Min Unbalanced chemistry integration time: "<<maxUnbalancedCPUtime<<"/"<<minUnbalancedCPUtime<<endl;
        Info<<"Max/Min Balanced execution time: "<<maxbalancedCPUtime<<"/"<<minbalancedCPUtime<<endl;
    }

    return deltaTMin;
}