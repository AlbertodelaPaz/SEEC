function [Hback,Hfoward]=Stokes(Trajectory,e,h,ToCompare)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Stokes Shift calculation
%	INPUT:
%	Trajectory      -Matrix of evolve sequences (M x L), every row
%                    corresponds to a different evolutionary step.
%		e           -Coupling matrix (qLxqL) formed from L x L submatrices
%			         of size q x q with the pairwise couplings for each pair 
%			         of sites (L submatrices in diagonal are irrelevant).
%		h           -fields matrix (q x L) with the single site contributions
%			         to the Hamiltonian.
%		ToCompare 	-Index of generations to be compare to the first one in
%                    the Trajectory variable.
%	OUTPUT:
%       Hback       -Vector of the hamiltonian back change.
%       Hfoward     -Vector of the hamiltonian foward change.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tries=length(ToCompare);
    L=size(Trajectory,2);
    Hback=zeros(L,tries);
    Hfoward=zeros(L,tries);
    for i=1:tries
        TestSeq=Trajectory(ToCompare(i),:);
        Native=Trajectory(1,:);
        Native2Test=repmat(Native,L,1);
        Test2Native=repmat(TestSeq,L,1);
        
        for j=1:L
            Native2Test(j,j)=TestSeq(1,j);
            Test2Native(j,j)=Native(1,j);
        end
        
        Hback(:,i)=Generalhamiltonian(Test2Native,e,h,2,1);
        Hback(:,i)=Hback(:,i)-Generalhamiltonian(TestSeq,e,h,2,1);
        Hfoward(:,i)=Generalhamiltonian(Native2Test,e,h,2,1);
        Hfoward(:,i)=Hfoward(:,i)-Generalhamiltonian(Native,e,h,2,1);
        
    end
end
