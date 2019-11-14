function [expS,aver]=effAlphabet(Trajectory,e,h)

[M,L]=size(Trajectory);
expS=zeros(M,L);
aver=zeros(M,1);


for i=1:M
    for j=1:L
        P=siteprobdistribution(Trajectory(i,:),h,e,j,1);
        expS(i,j)=exp(-sum(P.*log(P)));
    end
    aver(i)=mean(expS(i,:));
end

end