function Real_RFolder(path_data,M)

Nativ_index=[8940 13549 7670 58282 2976 65068 26059 78096 1 92771 77814];


if ~contains(pwd, path_data)
    cd(strcat(path_data, '/'));
end

addpath('/net/ab/cb/69/jxd162630/env.d/RestFamily')

Directory=dir('Alignment*.txt');
NoFamilies=size(Directory,1);

for i=1:NoFamilies
    %% Create directory
    Name_directory=strcat('Family', pad(num2str(i),3,'left','0'));
    Current_directory=[pwd '/' Name_directory '/'];
    
    align=0;
    load([Current_directory 'Parameters.mat']);
    fprintf(1, ['Family_', pad(num2str(i), 3,'left','0'), ' being analized. \n'] )
    
    
    
    
    %% Single evolution
    % (couplings)
    fprintf ( 1, ' Single Evolution: \n');
    fprintf ( 1, '  Calculating single trajectory (couplings). \n');
    Native=align(Nativ_index(i),:);
    [EvTraj,~,~,~,Susttime,R,Susttimeline,~]=Probevolution(Native,e,h,M);

    % (no couplings)
    fprintf ( 1, '  Calculating single trajectory (no couplings). \n');
    [EvTraj_null,~,~,~,Susttime_null,R_null,Susttimeline_null,~]=Probevolution(Native,zeros(size(e)),h,M);

    save([Current_directory 'R.mat'],'R');
    save([Current_directory 'R_null.mat'],'R_null');
    
    
    %% Overdispersion and Hamiltonian Plots
    % Single
    
    fprintf ( 1, ' Creating Figures \n');
    fprintf ( 1, '  Creating Overdispersion Trajectories\n');
    
    
    f3=figure(3);
    set(f3,'Position',[100 100 900 600])
    plot(1:50:M,R,'r','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Dispersion (R)')
    saveas(f3,[Current_directory '/Figures/R'],'png');
    saveas(f3,[Current_directory '/Figures/R']);
    
    
    f5=figure(5);
    set(f5,'Position',[100 100 900 600])
    plot(1:50:M,R_null,'r','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Dispersion (R)')
    saveas(f5,[Current_directory '/Figures/R_null'],'png');
    saveas(f5,[Current_directory '/Figures/R_null']);
    
    close(f3);
    close(f5);
  
    
    fprintf ( 1, ' Finished. \n');
    
    
    %% Average evolution
    fprintf ( 1, ' Average Evolution: \n');
    fprintf ( 1, '  Calculating average trajectories. \n');
    
    R=zeros(M/50,100);
    R_null=zeros(M/50,100);
    
    for j=1:100
    [EvTraj,~,~,~,~,R(:,j),~,~]=Probevolution(Native,e,h,M);
    [EvTraj_null,~,~,~,~,R_null(:,j),~,~]=Probevolution(Native,zeros(size(e)),h,M);
    end
    
    meanR=mean(R,2);
    sigmaR=std(R,0,2);
    meanR_null=mean(R_null,2);
    sigmaR_null=std(R_null,0,2);
    
    save([Current_directory 'Average_R.mat'],'meanR','sigmaR');
    save([Current_directory 'Average_R_null.mat'],'meanR_null','sigmaR_null');
    
    fprintf ( 1, '   Finished. \n');
    
    %% Overdispersion and Hamiltonian Average Plots
    % Single
    
    fprintf ( 1, '  Creating Average Overdispersion\n');

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
    saveas(f4,[Current_directory '/Figures/Average_R'],'png');
    saveas(f4,[Current_directory '/Figures/Average_R']);
    
    f5=figure(5);
    set(f5,'Position',[100 100 900 600])
    plot(1:50:M,meanR_null,'r','LineWidth',4)
    hold on
    plot(1:50:M,meanR_null+sigmaR_null,'b:','LineWidth',4)
    plot(1:50:M,meanR_null-sigmaR_null,'b:','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Dispersion (R)')
    saveas(f5,[Current_directory '/Figures/Average_R_null'],'png');
    saveas(f5,[Current_directory '/Figures/Average_R_null']);
    
    close(f4);
    close(f5);
    
    fprintf ( 1, '   Finished. \n');
    
    
    clearvars -except path_data Nativ_index M Directory NoFamilies i
    clc
end
end


% for i=1:11
% rmdir(strcat('Family', paset(gca,'FontSize',16)d(num2str(i),3,'left','0')));
% end
