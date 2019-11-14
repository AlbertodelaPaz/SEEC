% lecture les spectres pompe-sonde
% npassage permet de choisir le nombre de passages que l'on souhaite lire
% clear all;
path_data='\\CASSIS\Groupes\resphy\manips\datas\2018juin\fap3\pompesonde2ps\2018-06-22-17h35m44s\';
npassage=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = fopen([path_data 'Axe fréquentiel'],'r');
taille_sigma=fread(f1,1,'int32');
sigma = fread(f1,taille_sigma,'double')';
fclose(f1);
 
f2=fopen ([path_data 'T0'],'r');
taille_T0=fread(f2,1,'int32');
T0=fread(f2,taille_T0,'double');
fclose(f2);

f4=fopen ([path_data 'Bacground CCD'],'r');
taille_background=fread(f4,1,'int32');
Background_CCD =fread(f4,taille_background,'double');
fclose(f4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% lecture de tous les passages individuels

for ipassage=1:1:npassage
    if ipassage<10
        passage=strcat('passage00',int2str(ipassage));
    OFF_passage=strcat('OFF_passage00',int2str(ipassage));
    sonde=strcat('sonde00',int2str(ipassage));
    else
        passage=strcat('passage0',int2str(ipassage));
    OFF_passage=strcat('OFF_passage0',int2str(ipassage));
    sonde=strcat('sonde0',int2str(ipassage));
    end
    
    g1=fopen ([path_data passage],'r');
    g2=fopen ([path_data OFF_passage],'r');
    g3=fopen ([path_data sonde],'r');
    taille_spectre_sonde=fread(g3,1,'int32');
    
    pompe_sonde(:,:,ipassage)=fread(g1,[taille_sigma,taille_T0],'double');
    OFF_pompe_sonde(:,:,ipassage)=fread(g2,[taille_sigma,taille_T0],'double');
    Sonde(:,ipassage)=fread(g3,taille_sigma,'double'); % Sonde est moyenné sur tous les délais
    
    fclose(g1); 
    fclose(g2);
    fclose(g3);
end