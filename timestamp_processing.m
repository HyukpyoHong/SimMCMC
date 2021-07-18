%% For simulation data
cd('D:\OneDrive - kaist.ac.kr\Research\ResearchMaterial_HHP\TimeDelayEstimation\CodeForSobolev\Data_from_simulation\')
filename_list = dir('Priors_*');

filename_list.name

length(filename_list)

savelist = table('Size', [length(filename_list),1], 'VariableTypes', ["string"]);
savelist = cell(length(filename_list), 1);


for jj = 1:length(filename_list)
    savelist{jj} = filename_list(jj).name;
end
savetable = cell2table(savelist);

writetable(savetable, "filenamelist.csv", 'WriteRowNames', false, 'WriteVariableNames', false)


%% For real data
cd('D:\OneDrive - kaist.ac.kr\Research\ResearchMaterial_HHP\TimeDelayEstimation\CodeForSobolev\Data_from_simulation\KM_fix\')
filename_list = dir('Priors_*');

filename_list.name

length(filename_list)

savelist = table('Size', [length(filename_list),1], 'VariableTypes', ["string"]);
savelist = cell(length(filename_list), 1);


for jj = 1:length(filename_list)
    savelist{jj} = filename_list(jj).name;
end
savetable = cell2table(savelist);

writetable(savetable, "filenamelist_KM_fix.csv", 'WriteRowNames', false, 'WriteVariableNames', false)


%% For real data delay X pri
clear;clc;
% cd('D:\OneDrive - kaist.ac.kr\Research\ResearchMaterial_HHP\TimeDelayEstimation\CodeForSobolev\realdata_prior\')
cd('/Users/hyukpyohong/Dropbox/Twostep_delay/realdata_delayXshare_20210521/')
filename_list = dir('Priors_*');

filename_list.name

length(filename_list)

savelist = table('Size', [length(filename_list),1], 'VariableTypes', ["string"]);
savelist = cell(length(filename_list), 1);


for jj = 1:length(filename_list)
    savelist{jj} = filename_list(jj).name;
end
savetable = cell2table(savelist);

writetable(savetable, "filenamelist_realdata_delayXshare.csv", 'WriteRowNames', false, 'WriteVariableNames', false)


%% For real data delay Y fix 
clear;clc;
% cd('D:\OneDrive - kaist.ac.kr\Research\ResearchMaterial_HHP\TimeDelayEstimation\CodeForSobolev\realdata_prior\')
cd('D:/Dropbox/Dropbox/Twostep_delay/realdata_basic')
filename_list = dir('Theta_delY*');

filename_list.name

length(filename_list)

savelist = table('Size', [length(filename_list),1], 'VariableTypes', ["string"]);
savelist = cell(length(filename_list), 1);


for jj = 1:length(filename_list)
    savelist{jj} = filename_list(jj).name;
end
savetable = cell2table(savelist);

writetable(savetable, "filenamelist_realdata_delayYfix.csv", 'WriteRowNames', false, 'WriteVariableNames', false)




