function [P]=siteprobdistribution(sequence,h,e,ms,type)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Site probability distribution for evolution scripts
%	INPUT:
%		S          - Complete sequence (1 x L) with an alphabet from 1 to q.
%		e          - Coupling matrix (qLxqL) formed from L x L submatrices
%			         of size q x q with the pairwise couplings for each pair 
%			         of sites (L submatrices in diagonal are irrelevant).
%		h          - h fields matrix (q x L) with the single site contributions
%			         to the Hamiltonian.
%		ms         - site to generate local probability distribution.
%       type       - 1: Probability distribution
%                    2: Cumulative probability distribution.
%
%	OUTPUTS:
%		P          - Site probability distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [q,N]=size(h);
    P=zeros(q,1);
    
    for s=1:(ms-1)
        P=P+e((s-1)*q+sequence(s),(ms-1)*q+(1:q))';
    end
    
    P=P+h(:,ms);
    
    for s=(ms+1):N
        P=P+e((ms-1)*q+(1:q),(s-1)*q+sequence(s));
    end
    P=exp(P);
    P=P/sum(P);
    
    if (type==2)
       P=cumsum(P);
    end
    
end
