function H=Generalhamiltonian(Sequences,couplings,localfields,Htype,N1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	General Hamiltonian H
%	INPUT:
%		Sequences   - Translated alignment according to the
%			      return_alignment function of the DCA
%			      framework.
%		couplings   - Corresponding coupling matrix as returned
%			      from the DCAparameters function of the DCA
%			      framework.
%		localfields - Corresponding h fields matrix as returned
%			      from the DCAparameters function of the DCA
%			      framework.
%		Htype 	    - Hamiltonian type:
%                      Type 1: Sums couplings over
%			                   pairs across two species while
%                      Type 2: Evaluates the complete Potts
%                              Hamiltonian
%		N1 	    - Lenght of the first species. This number is
%			      ignored in the Htype 2 case.
%	OUTPUT:
%		H 	    - Vector of the Hamiltonians of each sequences.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    q=size(localfields,1);
    [M, N]=size(Sequences);
    H=zeros(M,1);
    
    %Local Fields
    for seq=1:M
        for res=1:N
            H(seq)=localfields((Sequences(seq,res)),res)+H(seq);
        end
    end
   
    %Couplings
    switch Htype
        case 1
            %First residue in the first domain
            %Compute across 2 different species
              for seq=1:M
                  for res=1:N1
                      for pair=(N1+1):N
                          iindex=q*(res-1)+Sequences(seq,res);
                          jindex=q*(pair-1)+Sequences(seq,pair);
                          H(seq)=couplings(iindex,jindex)+H(seq);
                      end
                  end
              end
        case 2
            %Complete Hamiltonian
            for seq=1:M
                  for res=1:N
                      if (res<N)
                          for pair=(res+1):N
                              iindex=q*(res-1)+Sequences(seq,res);
                              jindex=q*(pair-1)+Sequences(seq,pair);
                              H(seq)=couplings(iindex,jindex)+H(seq);
                          end
                      end
                  end
            end
    end
    H=-H;
end
