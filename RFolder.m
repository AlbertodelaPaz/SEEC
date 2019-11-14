function RFolder(path_data,M)

%Nativ_index=[8940 13549 7670 58282 2976 65068 94035 78096 1 92771 77814];


if ~contains(pwd, path_data)
    cd(strcat(path_data, '/'));
end

addpath('/net/uu/nm/bi/jxd162630/RestFamily')

Directory=dir('Alignment*.txt');
NoFamilies=size(Directory,1);

meanDispersion=zeros(M/50,NoFamilies);

for i=1:NoFamilies
    %% Create directory
    Name_directory=strcat('Family', pad(num2str(i),3,'left','0'));
    Current_directory=[pwd '/' Name_directory '/'];
    
    %% Parameter estimation
    load([Current_directory 'Stats.mat']);
    load([Current_directory 'Average_Stats.mat']);
    
    f2=figure(2);
    set(f2,'Position',[100 100 900 600])
    plot(1:50:M,R,'r','LineWidth',4)
    hold on
    plot(1:50:M,meanR+sigmaR,'b:','LineWidth',4)
    plot(1:50:M,meanR-sigmaR,'b:','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Dispersion (R)')
    saveas(f2,strcat('RFigures/R_', pad(num2str(i),3,'left','0')),'png');
    saveas(f2,strcat('RFigures/R_', pad(num2str(i),3,'left','0')));
    
    close(f2);
    meanDispersion(:,i)=meanR;


    fprintf ( 1, ' Finished. \n');
    
    clearvars -except path_data M Directory NoFamilies i meanDispersion
    clc
end

    f2=figure(2);
    set(f2,'Position',[100 100 900 600])
    plot(1:50:M,meanDispersion','LineWidth',4)
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    ylabel('Dispersion (R)')
    hold off
    saveas(f2,'HamiltonianFigures/R_trajectories','png');
    saveas(f2,'HamiltonianFigures/R_trajectories');
    close(f2);
    
end


% for i=1:11
% rmdir(strcat('Family', paset(gca,'FontSize',16)d(num2str(i),3,'left','0')));
% end
