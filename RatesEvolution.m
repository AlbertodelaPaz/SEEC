function [Trajectory,sitecount,sustcount,Sustitutiontime,Non_Syn,Syn,mutatedsites]=RatesEvolution(S,e,h,M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Evolution Simulation to calculate evolutionary rates
%	INPUT:
%		S          - Starting sequence of lenght L for the evolution code
%			         translated to number array (1 x L) with an alphabet
%			         from 1 to q.
%		e          - Coupling matrix (qLxqL) formed from L x L submatrices
%			         of size q x q with the pairwise couplings for each pair 
%			         of sites (L submatrices in diagonal are irrelevant).
%		h          - h fields matrix (q x L) with the single site contributions
%			         to the Hamiltonian.
%
%		M          - Number of evolutionary steps (Use multiples of 50).
%
%	OUTPUTS:
%		Trajectory         -Matrix of evolve sequences (M x L), every row corresponding
%                           to a different evolutionary step.
%       sitecount          -Number of times a single site was sampled along the
%                           simulation.
%       sustcount          -Number of non-synonymous substitution every site
%                           undergoes along the simulation
%       symtimeline        -Number of times a single site was sampled along the
%                           simulation.
%       Sustitutiontime    -Number of non-synonymous substitution every site
%                           undergoes along the simulation.
%		Non_syn            -Matrix of step-stamps (M x L) marking the step and site
%                           where a non-synounymous substitution happens.
%       Syn                -Matrix of step-stamps (M x L) marking the step and site
%                           where a synounymous substitution happens.
%       mutatedsites       -Logical vector (1 x L) marking sites that underwent
%                           at least one non-synounymous substitution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N=size(S,2);
    Trajectory=zeros(M,N);
    Trajectory(1,:)=S;
    sitecount=zeros(1,N);
    sustcount=zeros(1,N);
    flag=zeros(M,1);
     %Timeline=poissrnd(5,1,M);
    Timeline=1:M;
    Non_Syn=zeros(M,N);
    Syn=zeros(M,N);
    for generation=2:M
        ms=randvar(cumsum(ones(1,N)/N),1);
        distP=siteprobdistribution(Trajectory(generation-1,:),h,e,ms,2);
        alpha=randvar(distP,1);
        Sprueba=Trajectory(generation-1,:);
        Sprueba(ms)=alpha;
        sitecount(ms)=sitecount(ms)+1;
        Trajectory(generation,:)=Sprueba;
        if (Trajectory(generation,ms)~=Trajectory(generation-1,ms))
            sustcount(ms)=sustcount(ms)+1;
            flag(generation)=1;
            Non_Syn(generation,ms)=Timeline(generation);
        else
            Syn(generation,ms)=Timeline(generation);
        end
    end
   Sustitutiontime=Timeline(flag==1);
   mutatedsites=find(sustcount~=0);
end

function X=randvar(P,n)
% returns a realization for a discreate random variable given its
% cumulative probability distribution.
    [~,X] = histc(rand(1,n),P);
    X=X+1;
end
