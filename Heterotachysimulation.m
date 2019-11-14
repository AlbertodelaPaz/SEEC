function [Irregularityperc,IrregularityperMutations,sites,Countmut]=Heterotachysimulation(S,e,h,M,nosimulations,Tolerance,type)
    %type 1 : Poisson interpolation
    %type 2 : Average rate

    L=size(S,2);
    rates=zeros(nosimulations,L);
    
    for i=1:nosimulations
        [~,~,~,~,~,~,Sustitutiontimeline,~,~]=Probevolution(S,e,h,M);
        rates(i,:)=Poissonrate(Sustitutiontimeline,type);
        Countmut(i)=length(find(Sustitutiontimeline~=0));
    end
    [Irregularityperc,sites]=anormality(rates, Tolerance);
    actual=sum(rates);
    actual=sum(actual(:)==0);
    IrregularityperMutations=0;	
	if actual~=0
    IrregularityperMutations=Irregularityperc*L/actual;
	end
end

function [porcen,q]=anormality(freqs, Tolerance)
    numseq=size(freqs,1);
    numres=size(freqs,2);
    aver=repmat(mean(freqs),numseq,1);
    desv=repmat(std(freqs),numseq,1)/sqrt(numseq);
    z=(freqs-aver)./desv;
    alpha=(1-Tolerance)/2;
    zth=norminv(alpha,0,1);
    hetero=abs(z)>abs(zth);
    hetero=sum(hetero);
    q=hetero>((1-Tolerance)*numseq);
    porcen=sum(q)/numres;
end
