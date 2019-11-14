function stokesFolder(path_data)

if ~contains(pwd, path_data)
    cd(strcat(path_data, '/'));
end

addpath('/net/uu/nm/bi/jxd162630/Evo_Scripts')

Directory=string(path_data);

NoFamilies=size(Directory,1);
for i=1:NoFamilies
    Name_directory=Directory;
    Current_directory=Name_directory;
    
    %% Parameter estimation
    load(strcat(Current_directory, '/Parameters_orig.mat'));
    load(strcat(Current_directory, '/Stats.mat'));

   
    %Stokes
    fprintf ( 1, '  Calculating Stokes Shift. \n');
    [Hback,Hfoward]=Stokes(EvTraj,e,h,[250 1500]);
    Identity250=100-100*pdist([EvTraj(1,:);EvTraj(250,:)],'hamming');
    Identity1500=100-100*pdist([EvTraj(1,:);EvTraj(1500,:)],'hamming');
    save(strcat(Current_directory, '/Stokes.mat'),'Hback','Hfoward','Identity250','Identity1500');
    fprintf ( 1, ' Finished. \n');

    
    fprintf ( 1, ' Creating Figures \n');
    
    f6=figure(6);
    set(f6,'Position',[100 100 500 1000])
    for j=1:2
    subplot(2,1,j)
         plot(-100:100,100:-1:-100,'k','LineWidth',4)
         hold on
         scatter(Hback(:,j),Hfoward(:,j),100,'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor','none');    
         set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',4)
         xlabel('H_{back}')
         ylabel('H_{foward}')
         hold off
         axis([-1 1 -1 1]*10)
    end
    
    saveas(f6,strcat(Current_directory, '/Figures/Stokes'),'epsc');
    saveas(f6,strcat(Current_directory, '/Figures/Stokes'));


    close(f6)
    
    fprintf ( 1, ' Finished. \n');
    
    clearvars -except path_data M Directory NoFamilies i
    clc
end
end
