function HamiltonianFolder(path_data,M)

%Nativ_index=[8940 13549 7670 58282 2976 65068 94035 78096 1 92771 77814];


if ~contains(pwd, path_data)
    cd(strcat(path_data, '/'));
end

addpath('/net/ab/cb/69/jxd162630/env.d/LongerSimulations')

Directory=dir('Alignment*.txt');
NoFamilies=size(Directory,1);

meanHamiltonian=zeros(M,NoFamilies);

for i=1:NoFamilies
    %% Create directory
    Name_directory=strcat('Family', pad(num2str(i),3,'left','0'));
    Current_directory=[pwd '/' Name_directory '/'];
    
    %% Parameter estimation
    load([Current_directory 'Stats.mat']);
    load([Current_directory 'Average_Stats.mat']);
    
    f2=figure(2);
    set(f2,'Position',[100 100 900 600])
    plot(1:M,trajectoryH,'r','LineWidth',4)
    hold on
    plot(1:M,meanH+sigmaH,'b:','LineWidth',4)
    plot(1:M,meanH-sigmaH,'b:','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Hamiltonian')
    saveas(f2,strcat('HamiltonianFigures/Hamiltonian_trajectory_', pad(num2str(i),3,'left','0')),'png');
    saveas(f2,strcat('HamiltonianFigures/Hamiltonian_trajectory_', pad(num2str(i),3,'left','0')));
    
    close(f2);
    meanHamiltonian(:,i)=meanH;
    

    
    

    fprintf ( 1, ' Finished. \n');
    
    clearvars -except path_data M Directory NoFamilies i meanHamiltonian
    clc
end

    f2=figure(2);
    set(f2,'Position',[100 100 900 600])
    plot(1:M,meanHamiltonian','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Hamiltonian')
    hold off
    saveas(f2,'HamiltonianFigures/Hamiltonian_trajectories','png');
    saveas(f2,'HamiltonianFigures/Hamiltonian_trajectories');
    close(f2);
    
end


% for i=1:11
% rmdir(strcat('Family', paset(gca,'FontSize',16)d(num2str(i),3,'left','0')));
% end
