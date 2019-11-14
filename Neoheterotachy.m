function [heterotachy_degree, TN,TS,r, sites, finalrates, overallsites]=Neoheterotachy(Native,e,h,M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	1000 Evolution Simulations to calculate rates and heterotachy degree
%	INPUT:
%		Native     - Starting sequence of lenght L for the evolution code
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
%		  r                  -Cell with empirical mutation rate for sites with
%                           sitecount > 1 and sustcount > 1 in each simulation0
%       TS                 -Recurrence intervals for non-synonymous substituion
%                           for each site
%       TN                 -Recurrence intervals for sampling each site
%       sites              -Cell containing the indeces for the sites with
%                           sitecount > 1 and sustcount > 1 for each iteration
%       finalrates         -Empirical mutation rate for the sites with 
%                           sitecount > 1 and sustcount > 1 in all iterations
%       overallsites       -Indeces for which sitecount > 1 and sustcount > 1 in
%                           all iterations
%       heterotachy_degree -Heterotachy degree for each site.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=size(h,2);
sites=cell(1000,1);
TN=cell(1000,1);
TS=cell(1000,1);
r=cell(1000,1);

   fprintf ( 1, '  Evolving sequences. \n');
   
parfor i=1:1000
    [~,sitecount,sustcount,~,Non_Syn,Syn,~]=RatesEvolution(Native,e,h,M);
    [sites{i},TN{i},TS{i},r{i}]=Rates(Syn,Non_Syn,sitecount,sustcount);
end

fprintf ( 1, '  Obtainig rates. \n');

overallsites=sites{1};
for i=2:1000
overallsites=[overallsites sites{i}];
end
overallsites=unique(overallsites);
actualnumber=size(overallsites,2);
finalrates=zeros(1000,actualnumber);
for i=1:1000
    [~,iindex,jindex]=intersect(sites{i},overallsites);
    finalrates(i,jindex)=r{i}(iindex');
end

   fprintf ( 1, '  Estimating heterotachy degree. \n');

heterotachy_degree=zeros(1,L);
for i=1:actualnumber
    heterotachy_degree(overallsites(i))=length(find(isoutlier(finalrates(finalrates(:,i)~=0,i))));
end
heterotachy_degree=heterotachy_degree/10;
end
