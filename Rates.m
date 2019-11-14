function [sites,TN,TS,r]=Rates(Syn,Non_Syn,sitecount,sustcount)
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