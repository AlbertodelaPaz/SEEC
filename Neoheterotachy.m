function [heterotachy_degree, TN,TS,r, sites, finalrates, overallsites]=Neoheterotachy(Native,e,h,M)
L=size(h,2);
sites=cell(1000,1);
TN=cell(1000,1);
TS=cell(1000,1);
r=cell(1000,1);

   fprintf ( 1, '  Evolving sequences. \n');

parfor i=1:1000
    [~,sitecount,sustcount,~,~,Non_Syn,Syn,~]=RatesEvolution(Native,e,h,M);
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
