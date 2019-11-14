function Heterotachy_figures(path_data)

%Nativ_index=[8940 13549 7670 58282 2976 65068 94035 78096 1 92771 77814];
if ~contains(pwd,path_data)
    cd(strcat(path_data, '/'));
end

addpath('/net/uu/nm/bi/jxd162630/Evo_Scripts')

Directory=string(path_data);

NoFamilies=size(Directory,1);

for i=1:NoFamilies
    %% Create directory
    Name_directory=Directory;
    Current_directory=Name_directory;

    %% Parameter estimation
    load(strcat(Current_directory, '/Rates_.mat'));
    
    f2=figure(i);
    set(f2,'Position',[100 100 900 120])
    L=length(heterotachy_degree);
    imagesc(1,1:L,heterotachy_degree,[0,8])
    cb=colormap('hot');
    colormap(flip(cb))
    set(gca,'FontSize',18,'TickLength',[0.00 0.00],'LineWidth',3)
    xlabel('Sites')
    saveas(f2,strcat(Current_directory,'/Heterotachy'),'epsc');
   saveas(f2,strcat(Current_directory,'/Heterotachy'));
    
    close(f2);
    
    fprintf ( 1, ' Finished. \n');
    
    clearvars -except path_data Directory i
    clc
    end
    
end


% for i=1:11
% rmdir(strcat('Family', paset(gca,'FontSize',16)d(num2str(i),3,'left','0')));
% end
