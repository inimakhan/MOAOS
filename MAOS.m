
%__________________________________________________________________ %
%                          Multi-Objective                          %
%        Multi-Objetective Atomic Orbital Search (MOAOS)            %
%                                                                   %
%                                                                   %
%                  Developed in MATLAB R2022a (MacOs)               %
%                                                                   %
%                      Author and programmer                        %
%                ---------------------------------                  %
%                      Nima Khodadadi (ʘ‿ʘ)                         %
%                             e-Mail                                %
%                ---------------------------------                  %
%                         nkhod002@fiu.edu                          %
%                                                                   %
%                            Homepage                               %
%                ---------------------------------                  %
%                    https://nimakhodadadi.com                      %
%                                                                   %
%                                                                   %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------- %





%% Function Details
function [Archive_costs]=MAOS(MaxFes,Archive_size,nPop,VarNumber,method,m)

disp('MAOS Working');

%% General Parameters of the Algorithm

LayerNumber = 10 ;     % Maximum number of Layers around nucleus
FotonRate = 0.1 ;     % Foton Rate for position determination of electrons

%% Counters

% N=n_obj;
if method==3

    TestProblem=sprintf('P%d',m);

    CostFunction= Ptest(TestProblem);

    xrange  = xboundaryP(TestProblem);
    VarNumber=max(size(xrange));

    VarMin=xrange(:,1)';
    VarMax=xrange(:,2)';

end



%% Initialize the positions of solution




Alpha=0.1;  % Grid Inflation Parameter
nGrid=30;   % Number of Grids per each Dimension
beta=4; %=4;    % Leader Selection Pressure Parameter
gamma=2;

pop=zeros(nPop,VarNumber);
Pop=CreateEmptyParticle(nPop);


for i=1:nPop
    Pop(i).Velocity=0;
    Pop(i).Position=zeros(1,VarNumber);


    for j=1:VarNumber
        Pop(i).Position(1,j)=unifrnd(VarMin(j),VarMax(j),1);
    end
    pop(i,:)=Pop(i).Position;
    Pop(i).Cost=CostFunction(Pop(i).Position');

    Pop(i).Best.Position=Pop(i).Position;
    Pop(i).Best.Cost=Pop(i).Cost;
end

Pop=DetermineDominations(Pop);

Archive=GetNonDominatedParticles(Pop);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,Alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end
%Initialize the positions of moths

%% Initialization

MeanPop=mean(pop);

%% Main Loop
for FEs= 1:MaxFes-1

    NorDispCal=[];
    % Creat Quantum Layers
    MaxLay = randi(LayerNumber);
    NorDispInptut = 1:1:MaxLay;
    mu = 0;
    sigma = MaxLay/6;
    pd = makedist('Normal','mu',mu,'sigma',sigma);
    NorDisp = pdf(pd,NorDispInptut);
    NorDispCal(1,:)=NorDisp;
    NorDispCal(2,:)=NorDispCal(1,:)./sum(NorDispCal(1,:));
    NorDispCal(3,:)=nPop*NorDispCal(2,:);
    NorDispCal(4,:)=round(NorDispCal(3,:));
    NorDispCal(5,:)=cumsum(NorDispCal(4,:));
    LayCol=[0 NorDispCal(5,:)];
    LayCol(LayCol>nPop)=nPop;
    % Search Loop
    for i=1:MaxLay
        Leader1=SelectLeader(Archive,beta);
        Leader2=SelectLeader(Archive,beta);

        if size(Archive,1)>1
            counter=0;
            for newi=1:size(Archive,1)
                if sum(Leader1.Position~=Archive(newi).Position)~=0
                    counter=counter+1;
                    rep2(counter,1)=Archive(newi);
                end
            end
            Leader2=SelectLeader(rep2,beta);
        end
        for k=LayCol(i)+1:LayCol(i+1)
            PopA(k).Position=Pop(k , :).Position;
            PopA(k).Cost=Pop(k).Cost;
            popA(k,:)=Pop(k).Position;
            costA(k,:)=Pop(k).Cost;
            Pop(k).Position=zeros(1,VarNumber);
            Pop(k).Cost=0;

        end
        Energy=mean(costA);
        Orbit=i;
        for j=1:size(popA,1)
            if rand>FotonRate
                if costA(j,1)>Energy
                    Ir=unifrnd(0,1,1,2);
                    Jr=unifrnd(0,1,1,VarNumber);
                    Xold=popA(j,:);
                    Xbest=Leader1.Position;
                    Xmean=MeanPop;
                    Pop(j).Position=Xold+(Jr.*(Ir(1)*Xbest-Ir(2)*Xmean)/Orbit);
                    Pop(j).Position = max(Pop(j).Position,VarMin);
                    Pop(j).Position = min(Pop(j).Position,VarMax);
                    Pop(j).Cost=CostFunction(Pop(j).Position');

                else
                    Ir=unifrnd(0,1,1,2);
                    Jr=unifrnd(0,1,1,VarNumber);
                    Xold=popA(j,:);
                    Xbest=Leader2.Position;
                    if size(popA,1)==1
                        Xmean=popA;
                    else
                        Xmean=mean(popA);
                    end
                    Pop(j).Position=Xold+(Jr.*(Ir(1)*Xbest-Ir(2)*Xmean));
                    Pop(j).Position = max(Pop(j).Position ,VarMin);
                    Pop(j).Position = min(Pop(j).Position,VarMax);
                    Pop(j).Cost=CostFunction(Pop(j).Position');

                end
            else
                Pop(j).Position=unifrnd(VarMin,VarMax,[1 VarNumber]);
                Pop(j).Cost=CostFunction(Pop(j).Position');

            end
        end


    end

    FEs=FEs+1;
    Pop=DetermineDominations(Pop);
    non_dominated_Pop=GetNonDominatedParticles(Pop);

    Archive=[Archive
        non_dominated_Pop];

    Archive=DetermineDominations(Archive);
    Archive=GetNonDominatedParticles(Archive);

    for i=1:numel(Archive)
        [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end

    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);

        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,Alpha);

    end






    % Display the iteration and best optimum obtained so far
    disp(['In iteration ' num2str(FEs) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
    save results

    % Results

    costs=GetCosts(Pop);
    Archive_costs=GetCosts(Archive);



end




