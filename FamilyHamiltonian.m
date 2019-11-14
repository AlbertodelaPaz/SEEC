function FamilyHamiltonian(inputfile,M)

[h,e,align]=DCAparameters(inputfile,1);
N=size(align,1);     %Members
L=size(align,2);      %Length
Hinitial=zeros(M,N);
Hfinal=zeros(M,N);
CurrentTrajectory=zeros(M,L);
sustcount=zeros(1,L);
Sustitutions=zeros(1,N);
Final_align=zeros(N,L);

        Hinitial=Generalhamiltonian(align,e,h,2,1);

    parfor seq=1:N
        [CurrentTrajectory,~,sustcount,~,~,~,~,~,~]=Probevolution(align(seq,:),e,h,M);
        Hfinal(:,seq)=Generalhamiltonian(CurrentTrajectory,e,h,2,1);
        Sustitutions(seq)=sum(sustcount);
        Final_align(seq,:)=CurrentTrajectory(end,:);
        if (mod(seq,5000)==0)
            fprintf ( 1, 'Number of evolved sequences= %f\n', seq );
        end
    end

    save('Workspace_Whole_Family.mat');
    save('Dynamic_Variables.mat','Hfinal','align','Final_align','e','h','Hinitial');
    
end
