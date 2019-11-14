function HeterotachyFolder(path_data,M)

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
    align=0;
    Native=0;
    load(strcat(Current_directory, '/Parameters_orig.mat'));
    
    %Heterotachy and Rates
    fprintf ( 1, strcat(' Family ', Name_directory, '. \n'));
    
    fprintf ( 1, '  Calculating rates and heterotachy. \n');
    
    [heterotachy_degree, TN,TS,r, sites, finalrates, overallsites]=Neoheterotachy(Native,e,h,M);
    
    save(strcat(Current_directory, '/Rates_.mat'),'heterotachy_degree', 'TN','TS','r', 'sites', 'finalrates', 'overallsites');
    fprintf ( 1, ' Finished. \n');

    
    clearvars -except path_data M Directory NoFamilies i
    clc
end
end


% for i=1:11
% rmdir(strcat('Family', paset(gca,'FontSize',16)d(num2str(i),3,'left','0')));
% end
