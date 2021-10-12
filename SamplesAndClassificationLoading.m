% dir_old = 'C:\Users\SavirLab\OneDrive - Technion\TalBenYakov\SingleCellRNAseq\Fertility Immune\Data\Ovary_Old_2_Oct_19\';
% dir_old = 'C:\Users\SavirLab\Technion\Yoni Savir - TalBenYakov\SingleCellRNAseq\Fertility Immune\Data\Ovary_Old_2_Oct_19\';
dir_old = pwd;

data_10X_old = full(mmread([dir_old '\Data\Ovary_Old_2_Oct_19\matrix.mtx']));

[feature_ids_old gene_names_old feature_types_old] = textread([dir_old '\Data\Ovary_Old_2_Oct_19\features.tsv'],'%s %s %s','delimiter', '\t');

barcodes_old = textread([dir_old '\Data\Ovary_Old_2_Oct_19\barcodes.tsv'],'%s','delimiter', '\t');

%% Loading 10X data - Young ovary 26 Nov 19 

% dir_young = 'C:\Users\SavirLab\OneDrive - Technion\TalBenYakov\SingleCellRNAseq\Fertility Immune\Data\Ovary_Young_26_Nov_19\filtered\';
% dir_young = 'C:\Users\SavirLab\Technion\Yoni Savir - TalBenYakov\SingleCellRNAseq\Fertility Immune\Data\Ovary_Young_26_Nov_19\filtered\';
dir_young = pwd;

data_10X_young = full(mmread([dir_young '\Data\Ovary_Young_26_Nov_19\filtered\matrix.mtx']));

[feature_ids_young gene_names_young feature_types_young] = textread([dir_young '\Data\Ovary_Young_26_Nov_19\filtered\features.tsv'],'%s %s %s','delimiter', '\t');

barcodes_young = textread([dir_young '\Data\Ovary_Young_26_Nov_19\filtered\barcodes.tsv'],'%s','delimiter', '\t');

%% Load SingleR Classification results 
[~,txt_old,raw_old] = xlsread('classification_old.xlsx');
[~,txt_young,raw_young] = xlsread('classification_young.xlsx');

[~,txt_young_fine,raw_young_fine] = xlsread('classification_young_fine_2.xlsx');
[~,txt_old_fine,raw_old_fine] = xlsread('classification_old_fine.xlsx');