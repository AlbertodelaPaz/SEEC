function [sites,TN,TS,r]=Rates(Syn,Non_Syn,sitecount,sustcount)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Inference of the individual site rates. Input generated with RatesEvolution.m
%	INPUT:
%       sitecount          -Number of times a single site was sampled along the
%                           simulation.
%       sustcount          -Number of non-synonymous substitution every site
%                           undergoes along the simulation
%		Non_syn            -Matrix of step-stamps (M x L) marking the step and site
%                           where a non-synounymous substitution happens.
%       Syn                -Matrix of step-stamps (M x L) marking the step and site
%                           where a synounymous substitution happens.
%	OUTPUTS:
%		r                  -Empirical mutation rate for sites with sitecount > 1
%                           and sustcount > 1
%       TS                 -Recurrence intervals for non-synonymous substituion
%                           for each site
%       TN                 -Recurrence intervals for sampling each site
%       sites              -Index for the sites with sitecount > 1 and sustcount > 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sites=find((sitecount>1)&(sustcount>1));
    Syn=Syn(:,sites);
    Non_Syn=Non_Syn(:,sites);
    TN=zeros(size(sites));
    TS=zeros(size(sites));
    for i=1:length(sites)
        TS(i)=mean(diff(Syn(Syn(:,i)~=0,i)));
        TN(i)=mean(diff(Non_Syn(Non_Syn(:,i)~=0,i)));
    end
    r=TS./TN;
end
