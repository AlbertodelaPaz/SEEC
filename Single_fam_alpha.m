function Single_fam_alpha(path_data)
    
if ~contains(pwd, path_data)
    cd(strcat(path_data, '/'));
end

addpath('/net/uu/nm/bi/jxd162630/Evo_Scripts')

Directory=string(path_data);

NoFamilies=size(Directory,1);

%HeterSites=[78,13,24;... %Response 4
% 104,14,19];% Peripla 11

    
for i=1:NoFamilies
    %% Create directory
    Name_directory=Directory;
    Current_directory=Name_directory;
    
    %% Parameter estimation
    load(strcat(Current_directory, '/Stats.mat'));
    load(strcat(Current_directory, '/Parameters_orig.mat'));
    load(strcat(Current_directory, '/Rates_.mat'));
    
    [~,ind]=min(heterotachy_degree);
    HeterSites=ind;
    [~,ind]=min(abs(heterotachy_degree-median(heterotachy_degree)));
    HeterSites=[HeterSites ind];
    [~,ind]=max(heterotachy_degree);
    HeterSites=[HeterSites ind];
    
    [expS,averS]=effAlphabet(EvTraj(1:3000,:),e,h);
    
    f2=figure(2);
    set(f2,'Position',[100 100 1200 600])
    plot(1:3000,expS(:,HeterSites(1,1)),'Color',[1 0.3 0],'LineWidth',4)
    hold on
    plot(1:3000,expS(:,HeterSites(1,2)),'Color',[1 0.5 0],'LineWidth',4)
    plot(1:3000,expS(:,HeterSites(1,3)),'Color',[1 0.7 0],'LineWidth',4)
    hold off
    xlim([0 3000])
    set(gca,'FontSize',18,'TickLength',[0.03 0.03],'LineWidth',4)
    xlabel('Evolutionary Steps')
    legend('Low','Medium','High','Location','northeast')
    ylabel('Effective Alphabet')
    
    save(strcat(Current_directory, '/Alphabet_var.mat'),'expS','averS','HeterSites')
    saveas(f2,'Eff_Alph_bm','epsc');
    saveas(f2,'Eff_Alph_bm');
    
    
    close(f2);
    

    fprintf ( 1, ' Finished. \n');
    
    clearvars -except path_data M Directory NoFamilies i HeterSites
    clc
end


end
