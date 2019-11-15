function [expS,aver]=effAlphabet(Trajectory,e,h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Effective alphabet measure from the single-site probability distribution along trajectory
%	INPUT:
%		Trjectory  - Matrix of evolve sequences (M x L), every row corresponding
%                    to a different evolutionary step.
%		e          - Coupling matrix (qLxqL) formed from L x L submatrices
%			         of size q x q with the pairwise couplings for each pair 
%			         of sites (L submatrices in diagonal are irrelevant).
%		h          - h fields matrix (q x L) with the single site contributions
%			         to the Hamiltonian.
%
%	OUTPUTS:
%		expS       - Effective alphabet for every site along the trajectory
%		aver       - Average effective alphabet along the sites for every step.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
