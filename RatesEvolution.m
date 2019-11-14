function [Trajectory,sitecount,sustcount,Timeline,Sustitutiontime,Non_Syn,Syn,mutatedsites]=RatesEvolution(S,e,h,M)
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
