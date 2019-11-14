function GammaFolder(path_data)

%Nativ_index=[8940 13549 7670 58282 2976 65068 94035 78096 1 92771 77814];


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
    load(strcat(Current_directory, '/Rates_.mat'));
    
    gof=cell(1,1000);
    curve2=cell(1,1000);
    gr=cell(1,1000);
    for j=1:1000
    gr{j}=histogram(r{j}(r{j}<50),50,'Normalization','pdf');
    x=gr{j}.BinEdges(1:(end-1))+gr{j}.BinWidth/2;
    y=gr{j}.BinCounts/sum(gr{j}.BinCounts)/gr{j}.BinWidth;
    fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0],...
    'Upper',[Inf,Inf],...
    'StartPoint',[1 1]);
    ft = fittype('gampdf(x,a,b)','options',fo);
    [curve2{j},gof{j}] = fit(x',y',ft) ;
    adjusted(j)=gof{j}.adjrsquare;
    close(gcf)
    end
    
    [~,k] = max(adjusted);
        
     f2=figure(2);   
    set(f2,'Position',[100 100 900 600])
       hi=histogram(r{j}(r{j}<50),50,'Normalization','pdf');
        hold on
       x=hi.BinEdges(1:(end-1))+hi.BinWidth/2;
    y=hi.BinCounts/sum(hi.BinCounts)/hi.BinWidth;
     plot(x,gampdf(x,curve2{k}.a,curve2{k}.b),'LineWidth',3)
    set(gca,'FontSize',18,'TickLength',[0.02 0.02],'LineWidth',3)
    xlabel('$$\frac{T_S}{T_N}$$','Interpreter','latex')
    ylabel('Relative Frequency')
    legend({['Adjusted R-square: ' num2str(adjusted(k))]},'FontSize',16)
    hold off
    saveas(f2,strcat(Current_directory,'/Figures/GammaDist'),'epsc');
    saveas(f2,strcat(Current_directory,'/Figures/GammaDist'));
    close(f2);

    save(strcat(Current_directory,'/GammaDist.mat'))
    
    clearvars -except path_data M Directory NoFamilies i
    clc
end
    
end
