function siterate=Poissonrate(sustitutiontimeline,type)
    %type 1 : Poisson interpolation
    %type 2 : Average rate

    L=size(sustitutiontimeline,2);
    siterate=zeros(1,L);
    
    switch type
        case 1
    
    for i=1:L
        r=sustitutiontimeline(:,i);
        r=r(r~=0);
        if size(r,1)~=0
                r=diff(r);
            siterate(i)=poissfit(r);
        end
    end
    
        case 2
    for i=1:L
        r=sustitutiontimeline(:,i);
        r=r(r~=0);
        if size(r,1)~=0
                r=diff(r);
            siterate(i)=mean(r);
        end
    end

        case 3
    for i=1:L
        r=sustitutiontimeline(:,i);
        r=r(r~=0);
        if size(r,1)~=0
                r=diff(r);
            siterate(i)=mean(1./r);
        end
    end
end
end
