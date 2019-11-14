function evolveFolder(path_data,M)

if ~contains(pwd, path_data)
    cd(strcat(path_data, '/'));
end

addpath('/net/uu/nm/bi/jxd162630/Evo_Scripts')

Directory=string(path_data);

NoFamilies=size(Directory,1);

for i=1:NoFamilies
    Name_directory=Directory;
    Current_directory=Name_directory;
    mkdir(Name_directory);
    mkdir(strcat(Name_directory, '/Figures'))
    alignment=0;
    fprintf(1, strcat('Family ', Name_directory, ' being analized. \n'))
    
    %% Parameter estimation
    fprintf ( 1, '  Retriving DCA parameters: \n');
    Native=0;
    load(strcat(Directory(i), '/Parameters.mat'));
    align=alignment;
    Native=align(1,:);
    fprintf ( 1, '   Finished. \n');
    
    
    %% Hamiltonian
    fprintf ( 1, '  Calculating Family Hamiltonian\n');
    familyH=Generalhamiltonian(alignment,e,h,2,1);
    f1=figure(1);
    set(f1,'Position',[100 100 900 600])
    histogram(familyH,50,'Normalization','probability','FaceColor','g')
    set(gca,'FontSize',18,'LineWidth',4)
    xlabel('Hamiltonian')
    ylabel('Relative Frequency')
    saveas(f1,strcat(Current_directory, '/Figures/Family_Hamiltonian'),'png')
    saveas(f1,strcat(Current_directory, '/Figures/Family_Hamiltonian'))
    close(f1)
    fprintf ( 1, '   Finished. \n');
    save(strcat(Current_directory, '/Parameters_orig.mat'),'h','e','align','familyH','Native')
    
    %% Single evolution
    % (couplings)
    fprintf ( 1, ' Single Evolution: \n');
    fprintf ( 1, '  Calculating single trajectory (couplings). \n');
    [EvTraj,~,~,~,symtimeline,Susttime,R,Syn,Non_Syn,~]=Probevolution(Native,e,h,M);
    trajectoryH=Generalhamiltonian(EvTraj,e,h,2,1);
    
    % (no couplings)
    fprintf ( 1, '  Calculating single trajectory (no couplings). \n');
    [EvTraj_null,~,~,~,symtimeline_null,Susttime_null,R_null,Syn_null,Non_Syn_null,~]=Probevolution(Native,zeros(size(e)),h,M);
    trajectoryH_null=Generalhamiltonian(EvTraj_null,zeros(size(e)),h,2,1);

    save(strcat(Current_directory, '/Stats.mat'),'familyH','EvTraj','trajectoryH','R','symtimeline','Susttime','Syn','Non_Syn');
    save(strcat(Current_directory, '/Stats_null.mat'),'EvTraj_null','trajectoryH_null','R_null','symtimeline_null','Susttime_null','Syn_null','Syn_null');
    
    %Stokes
    fprintf ( 1, '  Calculating Stokes Shift. \n');
    [Hback,Hfoward]=Stokes(EvTraj,e,h,[500 1000 2000 5000]);
    [Hback_null,Hfoward_null]=Stokes(EvTraj_null,zeros(size(e)),h,[500 1000 2000 5000]);
    save(strcat(Current_directory, '/Stokes.mat'),'Hback','Hfoward');
    save(strcat(Current_directory, '/Stokes_null.mat'),'Hback_null','Hfoward_null');
    fprintf ( 1, ' Finished. \n');
    
    %% Overdispersion and Hamiltonian Plots
    % Single
    
    fprintf ( 1, ' Creating Figures \n');
    fprintf ( 1, '  Creating Overdispersion and Hamiltonian Trajectories\n');
    
    f2=figure(2);
    set(f2,'Position',[100 100 900 600])
    plot(1:M,trajectoryH,'r','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Hamiltonian')
    saveas(f2,strcat(Current_directory, '/Figures/Hamiltonian_trajectory'),'epsc');
    saveas(f2,strcat(Current_directory, '/Figures/Hamiltonian_trajectory'));
    
    f3=figure(3);
    set(f3,'Position',[100 100 900 600])
    plot(1:50:M,R,'r','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Dispersion (R)')
    saveas(f3,strcat(Current_directory, '/Figures/R'),'epsc');
    saveas(f3,strcat(Current_directory, '/Figures/R'));
    
    
    f4=figure(4);
    set(f4,'Position',[100 100 900 600])
    plot(1:M,trajectoryH_null,'r','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Hamiltonian')
    saveas(f4,strcat(Current_directory, '/Figures/Hamiltonian_trajectory_null'),'epsc');
    saveas(f4,strcat(Current_directory, '/Figures/Hamiltonian_trajectory_null'));
    
    f5=figure(5);
    set(f5,'Position',[100 100 900 600])
    plot(1:50:M,R_null,'r','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Dispersion (R)')
    saveas(f5,strcat(Current_directory, '/Figures/R_null'),'epsc');
    saveas(f5,strcat(Current_directory, '/Figures/R_null'));
    
    close(f2);
    close(f3);
    close(f4);
    close(f5);
    
    f6=figure(6);
    set(f6,'Position',[100 100 500 1000])
    for j=1:2
    subplot(2,1,j)
         scatter(Hback(:,j),Hfoward(:,j),100,'MarkerFaceColor',[0 .7 .7]);
         hold on
         plot(-100:100,100:-1:-100,'k','LineWidth',4)
         set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
         xlabel('H_{back}')
         ylabel('H_{foward}')
         hold off
         axis([-1 1 -1 1]*10)
    end
   
%     f7=figure(7);
%     set(f7,'Position',[100 100 500 1000])
%     for j=1:2
%     subplot(2,1,j)
%          scatter(Hback_null(:,j),Hfoward_null(:,j),'MarkerFaceColor',[0 .7 .7]);
%          hold on
%          plot(-10:10,10:-1:-10,'k','LineWidth',4)
%          set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
%          xlabel('H_{back}')
%          ylabel('H_{foward}')
%          hold off
%          axis([-1 1 -1 1]*20)
%     end
    
    saveas(f6,strcat(Current_directory, '/Figures/Stokes'),'epsc');
    saveas(f6,strcat(Current_directory, '/Figures/Stokes'));
%     saveas(f7,[Current_directory '/Figures/Stokes_null'],'png');
%     saveas(f7,[Current_directory '/Figures/Stokes_null']);

    close(f6)
%   close(f7)
    
    fprintf ( 1, ' Finished. \n');
    
    
    %% Average evolution
    fprintf ( 1, ' Average Evolution: \n');
    fprintf ( 1, '  Calculating average trajectories. \n');
    
    trajectoryH=zeros(M,100);
    R=zeros(M/50,100);
    trajectoryH_null=zeros(M,100);
    R_null=zeros(M/50,100);
    
    for j=1:100
    [EvTraj,~,~,~,~,~,R(:,j),~,~,~]=Probevolution(Native,e,h,M);
    trajectoryH(:,j)=Generalhamiltonian(EvTraj,e,h,2,1);
    [EvTraj_null,~,~,~,~,~,R_null(:,j),~,~,~]=Probevolution(Native,zeros(size(e)),h,M);
    trajectoryH_null(:,j)=Generalhamiltonian(EvTraj_null,zeros(size(e)),h,2,1);
    end
    
    meanH=mean(trajectoryH,2);
    sigmaH=std(trajectoryH,0,2);
    meanH_null=mean(trajectoryH_null,2);
    sigmaH_null=std(trajectoryH_null,0,2);
    
    meanR=mean(R,2);
    sigmaR=std(R,0,2);
    meanR_null=mean(R_null,2);
    sigmaR_null=std(R_null,0,2);
    
    save(strcat(Current_directory, '/Average_Stats.mat'),'meanH','sigmaH','meanR','sigmaR');
    save(strcat(Current_directory, '/Average_Stats_null.mat'),'meanH_null','sigmaH_null','meanR_null','sigmaR_null');
    
    fprintf ( 1, '   Finished. \n');
    
    %% Overdispersion and Hamiltonian Average Plots
    % Single
    
    fprintf ( 1, '  Creating Average Overdispersion and Hamiltonian Trajectories\n');
    
    f2=figure(2);
    set(f2,'Position',[100 100 900 600])
    plot(1:M,meanH,'r','LineWidth',4)
    hold on
    plot(1:M,meanH+sigmaH,'b:','LineWidth',4)
    plot(1:M,meanH-sigmaH,'b:','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Hamiltonian')
    hold off
    saveas(f2,strcat(Current_directory, '/Figures/Average_Hamiltonian_trajectory'),'epsc');
    saveas(f2,strcat(Current_directory, '/Figures/Average_Hamiltonian_trajectory'));
    
    f3=figure(3);
    set(f3,'Position',[100 100 900 600])
    plot(1:M,meanH_null,'r','LineWidth',4)
    hold on
    plot(1:M,meanH_null+sigmaH_null,'b:','LineWidth',4)
    plot(1:M,meanH_null-sigmaH_null,'b:','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Hamiltonian')
    saveas(f3,strcat(Current_directory, '/Figures/Average_Hamiltonian_trajectory_null'),'epsc');
    saveas(f3,strcat(Current_directory, '/Figures/Average_Hamiltonian_trajectory_null'));
    
    f4=figure(4);
    set(f4,'Position',[100 100 900 600])
    plot(1:50:M,meanR,'r','LineWidth',4)
    hold on
    plot(1:50:M,meanR+sigmaR,'b:','LineWidth',4)
    plot(1:50:M,meanR-sigmaR,'b:','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Dispersion (R)')
    hold off
    saveas(f4,strcat(Current_directory, '/Figures/Average_R'),'epsc');
    saveas(f4,strcat(Current_directory, '/Figures/Average_R'));
    
    f5=figure(5);
    set(f5,'Position',[100 100 900 600])
    plot(1:50:M,meanR_null,'r','LineWidth',4)
    hold on
    plot(1:50:M,meanR_null+sigmaR_null,'b:','LineWidth',4)
    plot(1:50:M,meanR_null-sigmaR_null,'b:','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Dispersion (R)')
    saveas(f5,strcat(Current_directory, '/Figures/Average_R_null'),'epsc');
    saveas(f5,strcat(Current_directory, '/Figures/Average_R_null'));
    
    close(f2);
    close(f3);
    close(f4);
    close(f5);
    
    fprintf ( 1, '   Finished. \n');
    
    %% Heterotachy Simulation
    %fprintf ( 1, '  Heterotachy Simulations. \n');
    %[Irregularityperc,IrregularityperMutations,sites,Countmut]=Heterotachysimulation(Native,e,h,M,100,.95,3);
    %[Irregularityperc_null,IrregularityperMutations_null,sites_null,Countmut_null]=Heterotachysimulation(Native,zeros(size(e)),h,M,100,.95,3);
    
    %save([Current_directory 'Heterotachy.mat'],'Irregularityperc','IrregularityperMutations','sites','Countmut');
    %save([Current_directory 'Heterotachy_null.mat'],'Irregularityperc_null','IrregularityperMutations_null','sites_null','Countmut_null');
    
    %fprintf ( 1, '   Finished. \n\n');
    
    clearvars -except path_data M Directory NoFamilies i
    clc
end
end


% for i=1:11
% rmdir(strcat('Family', paset(gca,'FontSize',16)d(num2str(i),3,'left','0')));
% end
