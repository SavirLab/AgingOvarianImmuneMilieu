%% Paper Graphs Script
close all
clc
clear

%% Loading Raw Data
SamplesAndClassificationLoading

%% Loading Suerat (R) parameters

[~,~,raw_cluster] = xlsread('All_tSNE_idnt_Jul21.csv'); % Cluster ID for each cell
[~,~,raw_tsne] = xlsread('All_tSNE_Jul21.csv'); % tSNE coordinates
[~,~,raw_age] = xlsread('All_tSNE_age_Jul21.csv'); % Sample (age) ID for each cell

all_tSNE = cell2mat (raw_tsne);
all_cluster= cell2mat(raw_cluster);
all_age = raw_age;

%% Excluding low quality cells by feature number and percentage of mitochindrial genes

data_both = {data_10X_young; data_10X_old};
mitoGenes = find(startsWith(gene_names_young,'mt-')); % find index of mitochondrial genes

for age=1:2 % for each sample
    nCounts{age} = sum(data_both{age}); % number of cells in sample
    
    for cell = 1:length(nCounts{age}) % for each cell in sample
        nFeatures{age}(cell) = length(find(data_both{age}(:,cell))); % how many features
        mitoFrction{age}(cell) = sum(data_both{age}(mitoGenes,cell))/sum(data_both{age}(:,cell)); % mitochondrial genes fraction
    end
    
    % We are using the following cutoffs: no less than 200 and no more than 2500 features & less than 10% mitochondrial genes
    data_reduced_ind{age} = find(nFeatures{age}>200 & nFeatures{age}<2500 & mitoFrction{age}<0.1); %which cells are ok
    data_reduced{age} = data_both{age}(:,data_reduced_ind{age}); % filtered data
end

%% Removing genes with no expression through both samples ("0 genes")

data_both_reduced = {data_reduced{1} data_reduced{2}};
data_both_joint = [data_reduced{1} data_reduced{2}];
gene_sum_joint = sum(data_both_joint,2); %which genes have any expression across all samples and cells 
expressed_genes_joint = find(gene_sum_joint); % index of expressed genes

% Removing "0 genes"
for age=1:2
    data_both_expressed{age} = data_both_reduced{age}(expressed_genes_joint,:);
end

data_both_expressed_joint = [data_both_expressed{1} data_both_expressed{2}];
%% Cluster annotation and Populations Setup based on Seurat results


color = [1 .4 .4; 1 .6 .2; .75 .75 .75;...
    .6 .6 0; .3 .6 0; .6 .3 0; 0 .6 .3; 0 .6 .6; 0 .8 .8;...
    .2 .6 1; .4 .7 1; .6 .6 1; 1 .6 1; 1 .2 .6]; % 0 .6 0; 


for cell = 1:length(all_age)
    if strcmp(all_age{cell},'young')
        age_num(cell) = 1;
    else
        age_num(cell) = 2;
    end
end

% excluding residual cells
all_cluster(all_cluster==15) = 0;
all_cluster(all_cluster==14) = 0;

%%%%%%%%%% For first-time analysis: %%%%%%%%%%%
% Setting all clusters except T cells (cluster 3)
% for age=1:2
% poptype_10X_R{1,age} = intersect(find(all_cluster==5), find(age_num==age))-(age-1)*sum(age_num==1); % Neutrophils
% poptype_num(1,age) = length(poptype_10X_R{1,age}); 
% 
% poptype_10X_R{2,age} = [intersect(find(all_cluster==7),find(age_num==age));...
%     intersect(find(all_cluster==10),find(age_num==age))]-(age-1)*sum(age_num==1); % Macrophages (7-M1?; 10-M2?)
% poptype_num(2,age) = length(poptype_10X_R{2,age});
% 
% poptype_10X_R{3,age} = [intersect(find(all_cluster==9),find(age_num==age));...
%     intersect(find(all_cluster==11),find(age_num==age))]-(age-1)*sum(age_num==1); % DC (11-cDC1; 9-cDC2)
% poptype_num(3,age) = length(poptype_10X_R{3,age}); 
% 
% poptype_10X_R{4,age} = intersect(find(all_cluster==6), find(age_num==age))-(age-1)*sum(age_num==1); % NK
% poptype_num(4,age) = length(poptype_10X_R{4,age}); 
% 
% poptype_10X_R{5,age} = intersect(find(all_cluster==1), find(age_num==age))-(age-1)*sum(age_num==1); % ILC1
% poptype_num(5,age) = length(poptype_10X_R{5,age}); 
% 
% poptype_10X_R{6,age} = intersect(find(all_cluster==4), find(age_num==age))-(age-1)*sum(age_num==1); % NKT
% poptype_num(6,age) = length(poptype_10X_R{6,age}); 
% 
% poptype_10X_R{9,age} = intersect(find(all_cluster==2), find(age_num==age))-(age-1)*sum(age_num==1); % DNT
% poptype_num(9,age) = length(poptype_10X_R{9,age}); 
% 
% poptype_10X_R{10,age} = intersect(find(all_cluster==8), find(age_num==age))-(age-1)*sum(age_num==1); % B
% poptype_num(10,age) = length(poptype_10X_R{10,age}); 
% 
% poptype_10X_R{11,age} = intersect(find(all_cluster==12), find(age_num==age))-(age-1)*sum(age_num==1); % B
% poptype_num(11,age) = length(poptype_10X_R{11,age}); 
% 
% poptype_10X_R{12,age} = intersect(find(all_cluster==13), find(age_num==age))-(age-1)*sum(age_num==1); % B
% poptype_num(12,age) = length(poptype_10X_R{12,age}); 
% 
% end
% 
% % Setting T cells and ILC2+3 clusters
% rest_pop = {'CD8','CD4'};
% tsne_age{1} = all_tSNE(age_num==1,:);
% tsne_age{2} = all_tSNE(age_num==2,:);
% 
% tsne_clst{1} = all_cluster(age_num==1,:);
% tsne_clst{2} = all_cluster(age_num==2,:);
% 
% for age = 1:2
% figure(0+age);
% scatter(tsne_age{age}(:,1),tsne_age{age}(:,2),10,'o','filled'); hold on
% 
% for T = 1:length(rest_pop)
% uiwait(msgbox(['Choose ' rest_pop{T} ' cells']))
%     h = impoly;
%     pol_vertices = h.getPosition;
%     ind_cut_type{T,age} = inpolygon(tsne_age{age}(:,1),tsne_age{age}(:,2),pol_vertices(:,1),pol_vertices(:,2));
%     
%     poptype_10X_R{6+T,age} = find(ind_cut_type{T,age}); % in data_map
%     scatter(tsne_age{age}(poptype_10X_R{6+T,age},1),tsne_age{age}(poptype_10X_R{6+T,age},2),5,'o', 'filled','MarkerFaceColor',color(T,:));    
%     poptype_num(6+T,age) = length(poptype_10X_R{6+T,age}); 
% end
% end
% 
% num_young = size(data_both_reduced{1},2);
% all_cluster(all_cluster==3) = 0;
% all_cluster([poptype_10X_R{7,1}; poptype_10X_R{7,2}+num_young]) = 3; % CD8
% all_cluster([poptype_10X_R{8,1}; poptype_10X_R{8,2}+num_young]) = 14; % CD4
% 
% all_cluster_temp = all_cluster;
% all_cluster_temp(all_cluster_temp==10) = 7;
% all_cluster_temp(all_cluster_temp==11) = 9;
% 
% save('R_10X_populations','poptype_10X_R','-v7.3');
% 
% cell_type_2clust = {'Neutrophils','Macrophsges','Dendritic','NK','ILC1','NKT','CD8 T cells','CD4 T cells','DNT cells', 'B cells','ILC2','ILC3'};

%%
%%%%%%%%%% For further analysis %%%%%%%%%%%
load('R_10X_populations.mat')

num_young = size(data_both_reduced{1},2);
all_cluster(all_cluster==3) = 0;
all_cluster([poptype_10X_R{7,1}; poptype_10X_R{7,2}+num_young]) = 3; % CD8
all_cluster([poptype_10X_R{8,1}; poptype_10X_R{8,2}+num_young]) = 14; % CD4

cell_type_2clust = {'Neutrophils','Macrophsges','Dendritic','NK','ILC1','NKT','CD8 T cells','CD4 T cells','DNT cells', 'B cells','ILC2','ILC3'};

all_cluster(all_cluster==10) = 7;

% all_cluster_temp = all_cluster;
% all_cluster_temp(all_cluster_temp==10) = 7;
% % all_cluster_temp(all_cluster_temp==11) = 9;

% ind_0 = find(all_cluster==0);
% data_both_joint_excel = data_both_joint;
% data_both_joint_excel(:,ind_0) = [];
% age_num_excel = age_num;
% age_num_excel(ind_0) = [];
% age_num_young = sum(age_num_excel==1);
% age_num_old = sum(age_num_excel==2);
% 
% all_cluster_reduced = all_cluster;
% all_cluster_reduced(ind_0) = [];
% 
% clst_order = [5 7 9 11 6 1 12 13 4 3 14 2 8];
% for i=1:13
%    ind_clst{i}=find(all_cluster_reduced==clst_order(i));
%    all_cluster_excel(ind_clst{i}) = i;
% end
% 
% ind_y = find(age_num_excel==1);
% ind_o = find(age_num_excel==2);
% data_excel = [all_cluster_excel' data_both_joint_excel'];
% data_excel_y = [all_cluster_excel(ind_y)' data_both_joint_excel(:,ind_y)'];
% data_excel_o = [all_cluster_excel(ind_o)' data_both_joint_excel(:,ind_o)'];
% writematrix(data_excel_o','data_old.csv')
% writematrix(data_excel_y','data_young.csv')
%%

%%%%%%%%%%%%%%%%%%
%   Figure 1B    %
%%%%%%%%%%%%%%%%%%

% ############# Should use the same color for both DC clusters??? #############

num_clst = max(all_cluster);

figure(3);
for clst =1:13%num_clst
scatter(all_tSNE(all_cluster==clst,1),all_tSNE(all_cluster==clst,2),5,color(clst,:),'o','filled'); hold on 
end

%% Violin plot for cluster identification

%%%%%%%%%%%%%%%%%%
%   Figure 1C    %
%%%%%%%%%%%%%%%%%%

%%%% Markers violin for cluster identification - 10X %%%%
ind_S100a8 = find(strcmp('S100a8',gene_names_young)); % Neutrophils
ind_Itgam = find(strcmp('Itgam',gene_names_young)); % Meyloids - Macs and DC
ind_Adgre1 = find(strcmp('Adgre1',gene_names_young)); % Macs 
ind_Tnf = find(strcmp('Tnf',gene_names_young)); % M1 
ind_Mrc1 = find(strcmp('Mrc1',gene_names_young)); % M2 
ind_Itgax = find(strcmp('Itgax',gene_names_young)); % DC
ind_Xcr1 = find(strcmp('Xcr1',gene_names_young)); % cDC1
ind_Cd209a = find(strcmp('Cd209a',gene_names_young)); % cDC2
ind_Klrb1c = find(strcmp('Klrb1c',gene_names_young)); % NK, ILC1 and NKT(?)
ind_Ccl5 = find(strcmp('Ccl5',gene_names_young)); % NK, ILC1, NKT T8  
ind_Itga1 = find(strcmp('Itga1',gene_names_young)); % ILC1 and NOT NK
ind_Gata3 = find(strcmp('Gata3',gene_names_young)); % ILC2 marker
ind_Il17rb = find(strcmp('Il17rb',gene_names_young)); % ILC2 marker
ind_Ly6c2 = find(strcmp('Ly6c2',gene_names_young)); %Mainly T8 or NKT
ind_Cd3e = find(strcmp('Cd3e',gene_names_young)); % Lymph
ind_Trbc2 = find(strcmp('Trbc2',gene_names_young)); % Lymph
ind_Cd8b1 = find(strcmp('Cd8b1',gene_names_young)); % Lymph
ind_Cd4 = find(strcmp('Cd4',gene_names_young));  % T4
ind_Cd28 = find(strcmp('Cd28',gene_names_young));  % T Lymphocytes
ind_Tmem176b = find(strcmp('Tmem176b',gene_names_young));  % DNT
ind_Il7r = find(strcmp('Il7r',gene_names_young)); % DNT + ILC2/3 > Ccr7+ T
ind_Tcrgc2 = find(strcmp('Tcrg-C2',gene_names_young)); % DNT and NKT and pre-T cells(?) 
ind_Il2ra = find(strcmp('Il2ra',gene_names_young)); % DNT, per-T, Tregs, ILC2/3 + T8/4  
ind_Cd19 = find(strcmp('Cd19',gene_names_young));  % B
ind_Ptprc = find(strcmp('Ptprc',gene_names_young));  % Pan immune


violin_markers = [ind_S100a8 ind_Itgam ind_Adgre1 ind_Itgax...
    ind_Klrb1c ind_Ccl5 ind_Itga1 ind_Gata3 ind_Il17rb ind_Ly6c2 ind_Cd3e ind_Trbc2 ind_Cd8b1 ind_Cd4 ind_Cd28...
    ind_Tmem176b ind_Il7r ind_Tcrgc2 ind_Il2ra ind_Cd19];
marker_name = {'S100a8','CD11b','F4/80','CD11c','NK1.1','Ccl5','CD49a',...
    'Gata3','Il17rb','Ly6c2','CD3e','TCRb-C2','CD8b1','CD4','CD28','Tmem176b','IL-7r','TCRg-C2','CD25','CD19'};

cluster_number = max(all_cluster);
num = 1:cluster_number;

for cnt_gene = 1:length(violin_markers)
    figure(4);
       subplot(20,1,cnt_gene); hold on
    
    gene_max = max(log(data_reduced{1}(violin_markers(cnt_gene),:)+1));
    gene_min = min(log(data_reduced{1}(violin_markers(cnt_gene),:)+1));
    k=1;
    for cnt_cluster = [5 7 9 11 6 1 12 13 4 3 14 2 8]
        
        data = data_reduced{1}(violin_markers(cnt_gene),all_cluster(age_num==1)==cnt_cluster);
        expressing = sum(data>0)/length(data);
        data_plot = ((log(data'+1)-gene_min)/(gene_max-gene_min))*expressing;
        violinplot_ys_lite(data_plot,num(k),0,'width',0.35,'ShowData',false,'ViolinColor',color(cnt_cluster,:),'MedianColor',color(cnt_cluster,:),'ShowMean',true);
        k=k+1;
        max_data_plot(cnt_gene,cnt_cluster) = max(data_plot);
    end
    
    top_value(cnt_gene) = round(max(max_data_plot(cnt_gene,:)),1);
    set(gca,'xtick',[1:13],'xticklabels',[])
    set(gca,'ytick',[0 top_value(cnt_gene)],'yticklabels',[]);
    xlim([0 13.5])
    ylabel(marker_name{cnt_gene})
end

%% SingleR classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Spplementary Figure 1     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_types_all = {'ILC (LIV.ILC1.DX5-)','ILC (ILC2)','ILC (ILC3.LTI.CD4-)','ILC (LPL.NCR+ILC3)',...
    'NK cells (NK.DAP10-)',...
    'NKT (NKT.4-)','NKT (NKT.44+NK1.1-)','NKT (NKT.44+NK1.1+)'...
    'Macrophages (MF.103-11B+24-)','Macrophages (MF.11C-11B+)','Macrophages (MFIO5.II+480INT)'...
    'Neutrophils (GN.ARTH)','Neutrophils (GN.Thio)',...
    'T cells (T.8EFF.OT1.D10LIS)','T cells (T.8MEM.OT1.D45.LISOVA)','T cells (T.CD4TESTCJ)','T cells (T.Tregs)','T cells (T.8EFF.TBET-.OT1LISOVA)',...
    'Tgd (Tgd.mat.VG1+VD6+)','Tgd (Tgd.mat.VG2+)',...
    'B cells (B.Fo)','B cells (B.FrF)','B cells (B.T3)','B cells (B1a)','B cells (B1A)',...
    'DC (DC.103+11B-)','DC (DC.103-11B+24+)','DC (DC.103-11B+F4-80LO.KD)','DC (DC.8+)','DC (DC.8-4-11B-)',...
    'Monocytes (MO.6C+II-)','Monocytes (MO.6C-II+)','Monocytes (MO.6C-IIINT)'};


txt_class_noHead_young_fine = txt_young_fine;
txt_class_noHead_old_fine = txt_old_fine(2:end);

for clst = 1:length(cell_types_all)
    type_young_all{clst} = find(contains(txt_class_noHead_young_fine(data_reduced_ind{1}),cell_types_all{clst}));
    type_old_all{clst} = find(contains(txt_class_noHead_old_fine(data_reduced_ind{2}),cell_types_all{clst}));
end

color_SingleR = [0.7 0 0; 1 0 0; 0.4 0 0; 0.4 0 0;...
    1 1 0;...
    1 0.7 0.4; 1 0.7 0.4; 1 0.7 0.4;...
    0.2 0.6 1; 0.2 0.6 1; 0.2 0.6 1;...
    0.2 1 0.2; 0.2 1 0.2;...
    .7 .7 .7; .7 .7 .7; 0 0 0; 0 0 0; .7 .7 .7;...
    1 0.2 1; 1 0.2 1;...
    0 1 1; 0 1 1; 0 1 1; 0 1 1; 0 1 1;...
    0.4 0 0.8; 0.4 0 0.8; 0.4 0 0.8; 0.4 0 0.8; 0.4 0 0.8;...
    0 0 1; 0 0 1; 0 0 1];

figure(5);
all_tSNE_Y = all_tSNE(age_num==1,:);
all_tSNE_O = all_tSNE(age_num==2,:);
scatter(all_tSNE(age_num==1,1),all_tSNE(age_num==1,2),5,'.','filled'); hold on
scatter(all_tSNE(age_num==2,1),all_tSNE(age_num==2,2),5,'.','filled');
for T = 1:length(type_young_all)
scatter(all_tSNE_Y(type_young_all{T},1),all_tSNE_Y(type_young_all{T},2),5,color_SingleR(T,:),'o','filled');
end
for T = 1:length(type_old_all)
scatter(all_tSNE_O(type_old_all{T},1),all_tSNE_O(type_old_all{T},2),5,color_SingleR(T,:),'o','filled'); 
end
% legend(['all' cell_types_all]);
title('SingleR FINE classifiction results on total 10X data - Young (no zeros) log(data+1) +PCA15')

%% CD8- CD4- in dofferent tissues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Spplementary Figure 3     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ############### for now only the young mouse id in the analysis #########

load('22_5.mat')
conditions_names = {'Ovary','Spleen','Peritoneum'};
conditions_names = categorical(conditions_names);
conditions_names = reordercats(conditions_names,{'Ovary','Spleen','Peritoneum'});

plotcells = cells_22_5_19_new([2,4,6],1:3);
% Columns 1:3 -> CD4, CD8, DN T cells (i.e CD3+ TCRb+)
% Rows 2,4,6 -> Young Ovaries, Spleen, Peritoneum

figure;
b = bar(conditions_names,plotcells,0.5,'stacked');
ylim([0 105])
b(1).FaceColor = [1 .4 .4];
b(2).FaceColor = [1 .6 .6];
b(3).FaceColor = [1 .8 .8];
legend('CD8','CD4','CD3+ TCRb+ CD4- CD8-')
box off


%%
%%%%%%%%%%%%%%%%%%
%   Figure 2A    %
%%%%%%%%%%%%%%%%%%

figure(6);
for age = 1:2
    for clst =1:cluster_number
        cells2plot = find(all_cluster==clst & age_num'==age);
        scatter3(all_tSNE(cells2plot,1),all_tSNE(cells2plot,2),age*ones(1,length(cells2plot)),...
            5,color(clst,:),'o','filled'); hold on
    end
end
set(gca,'xtick',[],'ytick',[],'ztick',[])

%% Joint tSNE by age
%%%%%%%%%%%%%%%%%%
%   Figure 2B    %
%%%%%%%%%%%%%%%%%%

figure(7);
scatter(all_tSNE(age_num==2,1),all_tSNE(age_num==2,2),10,[1 .6 .8],'o','filled','MarkerFaceAlpha',.5); hold on %[1 .6 .2]
scatter(all_tSNE(age_num==1,1),all_tSNE(age_num==1,2),10,[.4 .7 1],'o','filled','MarkerFaceAlpha',.5); 
title('tSNE on joint non-zero genes - log(Data+1)+PCA15')
legend({'Young','Old'});

%% Calculating Fractions + Bar plot
%%%%%%%%%%%%%%%%%%
%   Figure 2C    %
%%%%%%%%%%%%%%%%%%

for age=1:2
    for clst=1:size(poptype_10X_R,1)
        fractions(clst,age) = 100*(length(poptype_10X_R{clst,age})/sum(age_num==age));
        if isempty(fractions(clst,age))
            fractions(clst,age) =0;
        end
    end
    fractions(clst+1,age) = 100-sum(fractions(1:clst,age));
    poptype_num(13,age) = round((fractions(clst+1,age)/100)*sum(age_num==age)); 

end
cell_type_2clust{13} = 'Others';

age_group = {'Young','Old'};
figure(8); hold on
b=bar(fractions(:,:));
b(1).FaceColor = [.4 .7 1];
b(2).FaceColor = [1 .6 .8];
% legend(age_group);
set(gca,'xtick',1:13,'xticklabels',cell_type_2clust)
set(gca,'XTickLabelRotation',-30)
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels1 = string(poptype_num(:,1));
labels2 = string(poptype_num(:,2));

%% Figure 2; Fractions changes throughout age

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Figure 2D or supplementary 3    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Violin Plot - Facs (Macs+all CD3+ cells):
    % Go to Matlab folder and menually load all "cells" variables from the different
    % experiments (contains fractions of suppopulation).
load('25_7.mat');
load('19_12.mat');
load('16_1.mat');
load('9_8.mat');
    
cells_25_7_18(:,4) = [2 2 2 2 3 3 3];
cells_25_7_18(:,5) = [4 2 1 3 2 1 3];

cells_9_8_18(:,4) = [2 2 2 2 3 3 3 2];
cells_9_8_18(:,5) = [1 3 2 4 1 3 2 3];

cells_19_12_18(:,9) = [1 2 2 2 2 3 3 3 3];
cells_19_12_18(:,10) = [4 2 4 2 4 3 4 4 4];

for ii = 1:size(cells_19_12_18,1)
    cells_19_12_18(ii,7) = cells_19_12_18(ii,4) + cells_19_12_18(ii,5);
    cells_19_12_18(ii,11) = 100 - sum(cells_19_12_18(ii,[1,7]));
end
cells_19_12_4color = cells_19_12_18(:,[1 7 11 9:10]);

cells_16_1_19(:,9) = [1 1 1 2 2 3 3];
cells_16_1_19(:,10) = [1 4 2 3 4 4 4];
for ii = 1:size(cells_16_1_19,1)
    cells_16_1_19(ii,7) = cells_16_1_19(ii,4) + cells_16_1_19(ii,5);
    cells_16_1_19(ii,11) = 100 - sum(cells_16_1_19(ii,[1,7]));
end
cells_16_1_4color = cells_16_1_19(:,[1 7 11 9:10]);

cells_all = [cells_25_7_18; cells_9_8_18; cells_19_12_4color; cells_16_1_4color];

young_ind = 1;
adult_ind = 1;
old_ind = 1;
for ii = 1:size(cells_all,1)
    if cells_all(ii,4) == 1
        Macs_violin(young_ind,1) = cells_all(ii,1);
        Lymph_violin(young_ind,1) = cells_all(ii,2);
        young_ind = young_ind +1;
    elseif cells_all(ii,4) == 2
        Macs_violin(adult_ind,2) = cells_all(ii,1);
        Lymph_violin(adult_ind,2) = cells_all(ii,2);
        adult_ind = adult_ind +1;
    elseif cells_all(ii,4) == 3
        Macs_violin(old_ind,3) = cells_all(ii,1);
        Lymph_violin(old_ind,3) = cells_all(ii,2);
        old_ind = old_ind +1;
    end
end
Macs_violin(Macs_violin==0) = NaN;
Lymph_violin(Lymph_violin==0) = NaN;

cats = ({'Young','Adult','Old'});

for age = 1:3
    for stage = 1:4
        Macs_mean(stage,age) = nanmean(cells_all(cells_all(:,5)==stage & cells_all(:,4)==age,1));
        lymph_mean(stage,age) = nanmean(cells_all(cells_all(:,5)==stage & cells_all(:,4)==age,2));
    end
end

figure(9);
violins_macs = violinplot(Macs_violin, cats,'ViolinColor',[.3 .6 0],'Width',0.1,'ShowMean',true); hold on %[.4 .7 1]
violins_lymph = violinplot(Lymph_violin, cats, 'ViolinColor',[1 .7 .4],'Width',0.1,'ShowMean',true); %[1 .4 .4]
% run to add cycle stages:
plot(1:3,Macs_mean,'color',[.3 .6 0],'LineWidth',1.5)
plot(1:3,lymph_mean,'color',[1 .7 .4],'LineWidth',1.5)
%
% ylabel('Fraction from total CD45 cells (%)');

Young_sample = find(cells_all(:,4)==1);
Adult_sample = find(cells_all(:,4)==2);
Old_sample = find(cells_all(:,4)==3);

% Statistical analysis
% Macs
pval_FACS(1) = mattest(cells_all(Young_sample,1)',cells_all(Adult_sample,1)','VarType','unequal'); % Young vs. Adult
pval_FACS(2) = mattest(cells_all(Young_sample,1)',cells_all(Old_sample,1)','VarType','unequal'); % Young vs. Old
pval_FACS(3) = mattest(cells_all(Adult_sample,1)',cells_all(Old_sample,1)','VarType','unequal'); % Adult vs Old

%Lymph
pval_FACS(4) = mattest(cells_all(Young_sample,2)',cells_all(Adult_sample,2)','VarType','unequal'); % Young vs. Adult 
pval_FACS(5) = mattest(cells_all(Young_sample,2)',cells_all(Old_sample,2)','VarType','unequal'); % Young vs. Old
pval_FACS(6) = mattest(cells_all(Adult_sample,2)',cells_all(Old_sample,2)','VarType','unequal'); % Adult vs Old
%Macs vs Lymph
pval_FACS(7) = mattest(cells_all(Young_sample,1)',cells_all(Young_sample,2)','VarType','unequal'); % Young
pval_FACS(8) = mattest(cells_all(Adult_sample,1)',cells_all(Adult_sample,2)','VarType','unequal'); % Adut
pval_FACS(9) = mattest(cells_all(Old_sample,1)',cells_all(Old_sample,2)','VarType','unequal'); % Old


% Macs
[~,pval_FACS_ks(1)] = kstest2(cells_all(Young_sample,1)',cells_all(Adult_sample,1)'); % Young vs. Adult
[~,pval_FACS_ks(2)] = kstest2(cells_all(Young_sample,1)',cells_all(Old_sample,1)'); % Young vs. Old
[~,pval_FACS_ks(3)] = kstest2(cells_all(Adult_sample,1)',cells_all(Old_sample,1)'); % Adult vs Old

%Lymph
[~,pval_FACS_ks(4)] = kstest2(cells_all(Young_sample,2)',cells_all(Adult_sample,2)'); % Young vs. Adult 
[~,pval_FACS_ks(5)] = kstest2(cells_all(Young_sample,2)',cells_all(Old_sample,2)'); % Young vs. Old
[~,pval_FACS_ks(6)] = kstest2(cells_all(Adult_sample,2)',cells_all(Old_sample,2)'); % Adult vs Old
%Macs vs Lymph
[~,pval_FACS_ks(7)] = kstest2(cells_all(Young_sample,1)',cells_all(Young_sample,2)'); % Young
[~,pval_FACS_ks(8)] = kstest2(cells_all(Adult_sample,1)',cells_all(Adult_sample,2)'); % Adut
[~,pval_FACS_ks(9)] = kstest2(cells_all(Old_sample,1)',cells_all(Old_sample,2)'); % Old

%% Spider plot
%%%%%%%%%%%%%%%%%%
%   Figure 2E    %
%%%%%%%%%%%%%%%%%%

labels = {'Macrophages', 'CD11b^+ CD3^- cells', 'CD3^+ Lymphocytes', 'CD11b^- CD3^- cells'};

cells_spider = [cells_19_12_18; cells_16_1_19];
cells_spider_grouped(:,1) = cells_spider(:,1); 
cells_spider_grouped(:,2) = cells_spider(:,2)+cells_spider(:,3); 
cells_spider_grouped(:,3) = cells_spider(:,4)+cells_spider(:,5); 
cells_spider_grouped(:,4) = cells_spider(:,6)+cells_spider(:,8); 
cells_spider_grouped(:,5) = cells_spider(:,9); 

Young_spider = find(cells_spider_grouped(:,5)==1);
Adult_spider = find(cells_spider_grouped(:,5)==2);
Old_spider = find(cells_spider_grouped(:,5)==3);

spider_avg(1,:) = mean(cells_spider_grouped(Young_spider,:));
spider_std(1,:) = std(cells_spider_grouped(Young_spider,:));
spider_avg(2,:) = mean(cells_spider_grouped(Adult_spider,:));
spider_std(2,:) = std(cells_spider_grouped(Adult_spider,:));
spider_avg(3,:) = mean(cells_spider_grouped(Old_spider,:));
spider_std(3,:) = std(cells_spider_grouped(Old_spider,:));

% FACS data
figure(10); hold on
sp(1)=spider_plot(spider_avg(1,1:4), labels, 4, 1,'Marker', 'o','LineStyle', '-',...
    'LineWidth', mean(spider_std(1,1:4)),...
    'MarkerSize', 5, 'Color', [.4 .8 0],'MarkerFaceColor',[.4 .8 0]); %.4 .7 1
sp(2)=spider_plot(spider_avg(2,1:4), labels, 4, 1,'Marker', 'o','LineStyle', '-',...
    'LineWidth', mean(spider_std(2,1:4)),...
    'MarkerSize', 5,'Color', [1 .5 0],'MarkerFaceColor',[1 .5 0]); %.698 .4 1
sp(3)=spider_plot(spider_avg(3,1:4), labels, 4, 1,'Marker', 'o','LineStyle', '-',...
    'LineWidth', mean(spider_std(3,1:4)),...
    'MarkerSize', 5,'Color', [1 0 .5],'MarkerFaceColor',[1 0 .5]); %1 .6 .8 
sp(4)=spider_plot(spider_avg(1,1:4), labels, 4, 1,'Marker', 'o','LineStyle', '-',...
    'LineWidth', 2,...
    'MarkerSize', 5, 'Color', [.3 .6 0],'MarkerFaceColor',[.3 .6 0]);
sp(5)=spider_plot(spider_avg(2,1:4), labels, 4, 1,'Marker', 'o','LineStyle', '-',...
    'LineWidth', 2,...
    'MarkerSize', 5,'Color', [.8 .4 0],'MarkerFaceColor',[.8 .4 0]);
sp(6)=spider_plot(spider_avg(3,1:4), labels, 4, 1,'Marker', 'o','LineStyle', '-',...
    'LineWidth', 2,...
    'MarkerSize', 5,'Color', [.8 0 .4],'MarkerFaceColor',[.8 0 .4]);
for ii=1:3
    sp(ii).Color(4) = .8;
end

% 10X data
for age = 1:2
CD11b_noCD3(age) = fractions(1,age)+fractions(3,age)+fractions(4,age);
Cd3_noCD11b(age) = fractions(6,age)+fractions(7,age)+fractions(8,age)+fractions(9,age);
noCd3_noCD11b(age) = fractions(10,age)+fractions(11,age)+fractions(12,age)+fractions(5,age);
end

spider_10X = [fractions(2,1)  CD11b_noCD3(1) Cd3_noCD11b(1) noCd3_noCD11b(1);...
    fractions(2,2)  CD11b_noCD3(2) Cd3_noCD11b(2) noCd3_noCD11b(2)];

figure(11); hold on
sp10X(1)=spider_plot(spider_10X(1,1:4), labels, 4, 1,'Marker', 'o','LineStyle', '-',...
    'LineWidth', 4,...
    'MarkerSize', 5, 'Color', [.4 .8 0],'MarkerFaceColor',[.4 .8 0]); %.4 .7 1
sp10X(2)=spider_plot(spider_10X(2,1:4), labels, 4, 1,'Marker', 'o','LineStyle', '-',...
    'LineWidth', 4,...
    'MarkerSize', 5,'Color', [1 0 .5],'MarkerFaceColor',[1 0 .5]); %.698 .4 1

for ii=1:2
    sp10X(ii).Color(4) = .8;
end

%% Finding significantly changes genes

name_joint = gene_names_young(expressed_genes_joint); % Names of all expressed genes

% Calculating pVals including non-expressing cells
for clst = 1:size(poptype_10X_R,1)
    young{clst} = data_both_expressed{1}(:,poptype_10X_R{clst,1});
    old{clst} = data_both_expressed{2}(:,poptype_10X_R{clst,2});
    
    % significant DEGs (t-tesy)
    for gene = 1:length(expressed_genes_joint)
        [PValues{clst}(gene), TScores{clst}(gene)] = mattest(old{clst}(gene,:),young{clst}(gene,:),'VarType','unequal');
    end
    
    % Fraction diff
    frc_young{clst} = (sum(young{clst}>0,2)/size(young{clst},2))*100;
    frc_old{clst} = (sum(old{clst}>0,2)/size(old{clst},2))*100;
    FrcDiff{clst} = frc_old{clst} - frc_young{clst};
    
    % expression FC
    expression_level{clst,1} = mean(young{clst},2);
    nan_y = isnan(expression_level{clst,1});
    expression_level{clst,1}(nan_y) = 0;
    
    expression_level{clst,2} = mean(old{clst},2);
    nan_o = isnan(expression_level{clst,2});
    expression_level{clst,1}(nan_o) = 0;
    
    expression_FC{clst} = expression_level{clst,2}./expression_level{clst,1}; % mean(Old)/mean(Young)
    
    % Up or Downregulated?
    sig_up{clst} = intersect(find(PValues{clst}<=.05),find(log2(expression_FC{clst})>=1));
    sig_up_name{clst} = name_joint(sig_up{clst});
    sig_up_FC{clst} = expression_FC{clst}(sig_up{clst});
    sig_up_pval{clst}(:,1) = PValues{clst}(sig_up{clst});
    sig_down{clst} = intersect(find(PValues{clst}<=.05),find(log2(expression_FC{clst})<=-1));
    sig_down_name{clst} = name_joint(sig_down{clst});
    sig_down_FC{clst} = expression_FC{clst}(sig_down{clst});
    sig_down_pval{clst}(:,1) = PValues{clst}(sig_down{clst});
    
    % expression proportion change (chi-square)
    expression_fraction{clst,1} = 100*(sum(young{clst}>0,2)./size(young{clst},2));
    expression_fraction{clst,2} = 100*(sum(old{clst}>0,2)./size(old{clst},2));
    
    clst
end

%%
%%%%%%%%%%%%%%%%%%
%    Figure 3A    %
%%%%%%%%%%%%%%%%%%
% ############### Need to add DRGs per cluster menually ###############################

ind_Ifngr1 = find(strcmp('Ifngr1',name_joint)); % NT, Macs
ind_Anxa1 = find(strcmp('Anxa1',name_joint)); % NT
ind_Tgfbr1 = find(strcmp('Tgfbr1',name_joint)); % NT
ind_Cxcl2 = find(strcmp('Cxcl2',name_joint)); % NT,Macs,DC,NK,ILC1,NKT,CD8,CD4,DNT,B,ILC3
ind_Cxcr6 = find(strcmp('Cxcr6',name_joint)); % NT
ind_Il1a = find(strcmp('Il1a',name_joint)); % NT
ind_Il1b = find(strcmp('Il1b',name_joint)); % NT, ILC1, NKT
ind_Ptgs2 = find(strcmp('Ptgs2',name_joint)); % NT
ind_Cxcl10 = find(strcmp('Cxcl10',name_joint)); % NT, Macs, DNT
ind_Csf2rb = find(strcmp('Csf2rb',name_joint)); % NT
ind_Vegfa = find(strcmp('Vegfa',name_joint)); % NT, Macs

red_dots{1} = [ind_Ifngr1 ind_Anxa1 ind_Tgfbr1 ind_Cxcl2 ind_Cxcr6 ind_Il1a...
    ind_Il1b ind_Ptgs2 ind_Cxcl10 ind_Csf2rb ind_Vegfa];
red_dots_names{1} = {'Ifngr1','Anxa1','Tgfbr1','Cxcl2','Cxcr6','Il1a','Il1b',...
    'Ptgs2','Cxcl10','Csf2rb','Vegfa'};

ind_Il2ra = find(strcmp('Il2ra',name_joint)); % Macs
ind_F11r = find(strcmp('F11r',name_joint)); % Macs
ind_Gata3 = find(strcmp('Gata3',name_joint)); % Macs
ind_Mrc1 = find(strcmp('Mrc1',name_joint)); % Macs
ind_F10 = find(strcmp('F10',name_joint)); % Macs
ind_Igf1 = find(strcmp('Igf1',name_joint)); % Macs
ind_Hif1a = find(strcmp('Hif1a',name_joint)); % Macs
ind_Vcam1 = find(strcmp('Vcam1',name_joint)); % Macs
ind_Cd274 = find(strcmp('Cd274',name_joint)); % Macs
ind_Tnf = find(strcmp('Tnf',name_joint)); % Macs

red_dots{2} = [ind_Ifngr1 ind_Cxcl2 ind_Cxcl10 ind_Il2ra ind_F11r...
    ind_Gata3 ind_Mrc1 ind_F10 ind_Igf1 ind_Hif1a ind_Vcam1 ind_Cd274 ind_Tnf];
red_dots_names{2} = {'Ifngr1','Cxcl2','Cxcl10','Il2ra','F11r','Gata3','Mrc1',...
    'F10','Igf1','Hif1a','Vcam1','Cd274','Tnf'};

ind_Cd81 = find(strcmp('Cd81',name_joint)); % DC
ind_F2r = find(strcmp('F2r',name_joint)); % DC
ind_Ccr5 = find(strcmp('Ccr5',name_joint)); % DC, ILC1
ind_Cxcl16 = find(strcmp('Cxcl16',name_joint)); % DC
ind_Stat1 = find(strcmp('Stat1',name_joint)); % DC,NK,ILC1

red_dots{3} = [ind_Cd81 ind_Cxcl2 ind_F2r ind_Ccr5 ind_Cxcl16 ind_Stat1];
red_dots_names{3} = {'Cd81','Cxcl2','F2r','Ccr5','Cxcl16','Stat1'};

ind_Stat4 = find(strcmp('Stat4',name_joint)); % NK
ind_Icam1 = find(strcmp('Icam1',name_joint)); % NK
ind_Ccl3 = find(strcmp('Ccl3',name_joint)); % NK

red_dots{4} = [ind_Stat4 ind_Cxcl2 ind_Icam1 ind_Ccl3 ind_Stat1];
red_dots_names{4} = {'Stat4','Cxcl2','Icam1','Ccl3','Stat1'};

red_dots{5} = [ind_Cxcl2 ind_Il1b ind_Ccr5 ind_Stat1];
red_dots_names{5} = {'Cxcl2','Il1b','Ccr5','Stat1'};

ind_Ccl5 = find(strcmp('Ccl5',name_joint)); % NKT, ILC2
ind_Fcer1g = find(strcmp('Fcer1g',name_joint)); % NKT
ind_Mmp9 = find(strcmp('Mmp9',name_joint)); % NKT

red_dots{6} = [ind_Cxcl2 ind_Il1b ind_Ccl5 ind_Fcer1g ind_Mmp9];
red_dots_names{6} = {'Cxcl2','Il1b','Ccl5','Fcer1g','Mmp9'};

ind_Ly6c2 = find(strcmp('Ly6c2',name_joint)); % CD8, CD4

red_dots{7} = [ind_Cxcl2 ind_Ly6c2];
red_dots_names{7} = {'Cxcl2','Ly6c2'};

ind_Icam2 = find(strcmp('Icam2',name_joint)); % CD4, ILC2

red_dots{8} = [ind_Cxcl2 ind_Ly6c2 ind_Icam2];
red_dots_names{8} = {'Cxcl2','Ly6c2','Icam2'};

ind_Mmp16 = find(strcmp('Mmp16',name_joint)); % DNT
ind_Ccl4 = find(strcmp('Ccl4',name_joint)); % DNT,B

red_dots{9} = [ind_Cxcl2 ind_Mmp16 ind_Ccl4];
red_dots_names{9} = {'Cxcl2','Mmp16','Ccl4'};

ind_Iglc2 = find(strcmp('Iglc2',name_joint)); % B

red_dots{10} = [ind_Cxcl2 ind_Ccl4 ind_Iglc2];
red_dots_names{10} = {'Cxcl2','Ccl4','Iglc2'};

ind_Il10rb = find(strcmp('Il10rb',name_joint)); % ILC2
ind_Ccr2 = find(strcmp('Ccr2',name_joint)); % ILC2

red_dots{11} = [ind_Il10rb ind_Ccr2 ind_Icam2 ind_Ccl5];
red_dots_names{11} = {'Il10rb','Ccr2','Icam2','Ccl5'};

ind_Xcl1 = find(strcmp('Xcl1',name_joint)); % ILC3
ind_Il7r = find(strcmp('Il7r',name_joint)); % ILC3
ind_Ccr6 = find(strcmp('Ccr6',name_joint)); % ILC3
ind_Prdx2 = find(strcmp('Prdx2',name_joint)); % ILC3

red_dots{12} = [ind_Xcl1 ind_Il7r ind_Ccr6 ind_Prdx2 ind_Cxcl2];
red_dots_names{12} = {'Xcl1','Il7r','Ccr6','Prdx2','Cxcl2'};


for clst = 1:size(poptype_10X_R,1)
   x = log2(expression_FC{clst});
   y = -log10(PValues{clst});
    
   figure(11+clst);
   scatter(x,y,8,[.6 .6 .6],'o','filled'); hold on
   scatter(x(red_dots{clst}),y(red_dots{clst})',10,'ro','filled');
   text(x(red_dots{clst})+0.1,y(red_dots{clst})+0.1,red_dots_names{clst},'FontSize',13)
   line([1 1],[min(y)-0.5 max(y)+0.5],'LineStyle','--')
   line([-1 -1],[min(y)-0.5 max(y)+0.5],'LineStyle','--')
   line([min(x(x>-Inf))-0.5 max(x(x<Inf))+0.5],[0.585 0.585],'LineStyle','--')
   ylim([0 max(y)+0.5])
   xlim([-6 6])
   title(cell_type_2clust{clst})
end
    
%% 
%%%%%%%%%%%%%%%%%%
%    Figure 3B    %
%%%%%%%%%%%%%%%%%%

REVIGO_down_name = {'NT','Macs','DC','NK','ILC1','NKT','CD8','CD4','DNT','B','ILC2','ILC3'};
dir = 'C:\Users\SavirLab\Technion\Yoni Savir - TalBenYakov\Fertility Immune paper\REVIGO Jul21\';
k=2;
total_GO_list = {'GO term'};

for clst = 1:length(REVIGO_down_name)
    [~,txt_REVIGO_down{clst},raw_REVIGO_down{clst}] = xlsread([dir 'REVIGO_Down_' REVIGO_down_name{clst} '.csv']);
    for go = 2:size(raw_REVIGO_down{clst},1)
        if raw_REVIGO_down{clst}{go,11} == 0
            flg=0;
            for i=1:size(total_GO_list,1)
                if contains(raw_REVIGO_down{clst}{go,1},total_GO_list{i,1})
                    flg=1;
                end
            end
            if flg == 0   
                total_GO_list(k,1) = raw_REVIGO_down{clst}(go,1);
                total_GO_list(k,2) = raw_REVIGO_down{clst}(go,2);
                k=k+1;
            end
        end
    end
end

total_GO_list(1,3:length(REVIGO_down_name)+2) = REVIGO_down_name;
for term = 2:size(total_GO_list,1)
    n=1;
    for clst = 1:length(REVIGO_down_name)
        for i = 2:size(raw_REVIGO_down{clst},1)
            if contains(total_GO_list{term,1},raw_REVIGO_down{clst}{i,1})
                total_GO_list{term,clst+2} = 1;
                total_GO_list{term,12} = n;
                n=n+1;
            end
        end
    end
end

for r = 1:size(total_GO_list,1)
    for c = 1:size(total_GO_list,2)
      if isempty(total_GO_list{r,c})
          total_GO_list{r,c} = 0;
      end
    end
end


[~,txt_REVIGO_common_down,raw_REVIGO_common_down] = xlsread([dir 'REVIGO_Down_Common.csv']);
figure(24)
for i=2:size(raw_REVIGO_common_down,1)
    if isnumeric(raw_REVIGO_common_down{i,4})
        scatter(raw_REVIGO_common_down{i,4},raw_REVIGO_common_down{i,5},3^raw_REVIGO_common_down{i,6},raw_REVIGO_common_down{i,7},...
            'o','filled','MarkerEdgeColor',[0 0 0]); hold on
    end
end

for go = 2:size(total_GO_list,1)
    for clst  = 1:length(REVIGO_down_name)
        if total_GO_list{go,clst+2} == 1 && total_GO_list{go,end} == 1
            for row = 1:size(raw_REVIGO_down{clst},1)
                if cell2mat(strfind(raw_REVIGO_down{clst}(row,1),total_GO_list{go,1}))
                    raw_REVIGO_down{clst}{row,12} = 1;
                end
            end
        elseif total_GO_list{go,clst+2} == 1 && total_GO_list{go,end} > 1
            for row = 1:size(raw_REVIGO_down{clst},1)
                if cell2mat(strfind(raw_REVIGO_down{clst}(row,1),total_GO_list{go,1}))
                    raw_REVIGO_down{clst}{row,12} = 0;
                end
            end
        end
    end
end



%% Finding Cytokines and Chemokines 
chemokines_prefix = {'Ccl','Cxcl','Cx3cl','Xcl','Xcr','Ackr','Ccr','Cxcr','Cx3cr'};
ind=1;
for prefix = 1:length(chemokines_prefix)
chemokines_ind = find(startsWith(name_joint,chemokines_prefix{prefix}));
current_ind = ind:ind-1+length(chemokines_ind);
chemokines_idx(current_ind,1) = chemokines_ind;
chemokines(current_ind,1) = name_joint(chemokines_ind);
ind=ind+length(chemokines_ind);
end

cytokines_prefix = {'Amh','Bmp','Cd40lg','Cd70','Clcf1','Csf','Eda','Epo','Fas','Gdf','Ifn','Il','Inh',...
    'Lep','Lif','Lta','Ltb','Mstn','Ngf','Nodal','Osm','Prl','Tgf','Tnf','Tpo','Tslp','Acvr','Cd27','Cd4',...
    'Cntfr','Mpl','Relt'};
ind=1;
for prefix = 1:length(cytokines_prefix)
cytokines_ind = find(startsWith(name_joint,cytokines_prefix{prefix}));
current_ind = ind:ind-1+length(cytokines_ind);
cytokines_idx(current_ind,1) = cytokines_ind;
cytokines(current_ind,1) = name_joint(cytokines_ind);
ind=ind+length(cytokines_ind);
end

not_cyto = [21,23,24,26:30,41,45,57,59,66,74,84,86,87,92,101,113,119,121,124,...
    127:130,132,133,139,145,146,148,152,154,168,173,175,176,178,179,186,189,...
    199:203,205,207,209];

cytokines(not_cyto) = [];
cytokines_idx(not_cyto) = [];

%%
%%%%%%%%%%%%%%%%%%
%    Figure 4A   %
%%%%%%%%%%%%%%%%%%

clear l
clear r
clear G
clear c
[~,~,c] =xlsread('Chemokine network_KEGG_template.xlsx');

l = c(2:end,1);
r = c(2:end,7);
clear chemokines_map_Young
clear chemokines_map_Old
clear chemokines_map_h
clear ligands
clear receptors

for clst=1:12
    for gene = 1:length(chemokines)
        FC = mean(old{clst}(chemokines_idx(gene),:),2)/mean(young{clst}(chemokines_idx(gene),:),2);
       if  PValues{clst}(chemokines_idx(gene))<.05 && FC>=2 && FC<Inf
           if -log10(PValues{clst}(chemokines_idx(gene)))<=3
               chemokines_map(clst,gene) = 1;
           elseif -log10(PValues{clst}(chemokines_idx(gene)))>3 && -log10(PValues{clst}(chemokines_idx(gene)))<=6
               chemokines_map(clst,gene) = 2;
           elseif -log10(PValues{clst}(chemokines_idx(gene)))>6
               chemokines_map(clst,gene) = 3;
           end
       elseif PValues{clst}(chemokines_idx(gene))<.05 && FC<=.5 && FC>0
           if -log10(PValues{clst}(chemokines_idx(gene)))<=3
               chemokines_map(clst,gene) = -1;
           elseif -log10(PValues{clst}(chemokines_idx(gene)))>3 && -log10(PValues{clst}(chemokines_idx(gene)))<=6
               chemokines_map(clst,gene) = -2;
           elseif -log10(PValues{clst}(chemokines_idx(gene)))>6 && -log10(PValues{clst}(chemokines_idx(gene)))<=9
               chemokines_map(clst,gene) = -3;
           elseif -log10(PValues{clst}(chemokines_idx(gene)))>9
               chemokines_map(clst,gene) = -4;
           end
       else
           chemokines_map(clst,gene) = 0;
       end
    end
end

no_zero_clst_chmo = find(sum(chemokines_map~=0));
no_zero_chemokines = chemokines(no_zero_clst_chmo);

lig = 1;
rec = 1;
for clst=1:length(no_zero_chemokines)
   if sum(strcmp(no_zero_chemokines{clst},l))>0
    lig_flag_chmo(lig) = clst;
    lig=lig+1;
   elseif sum(strcmp(no_zero_chemokines{clst},r))>0
    rec_flag_chmo(rec) = clst;
    rec=rec+1;
   end
end
no_zero_chemokines_ordered = no_zero_chemokines([lig_flag_chmo 16 rec_flag_chmo]);
no_zero_clst_ordered_chmo = no_zero_clst_chmo([lig_flag_chmo 16 rec_flag_chmo]);


clr_map = [1 .2 .2; 1 .4 .4; 1 .6 .6; 1 .8 .8; 1 1 1; .8 1 .6]; %[1 .2 .2; 1 .4 .4; 1 .6 .6; 1 .8 .8; 1 1 1; .9 1 .8; .8 1 .6; .7 1 .4]%; .6 1 .2;];
figure(25);
heatmap(chemokines_map(:,no_zero_clst_ordered_chmo),'Colormap',clr_map,'ColorbarVisible','off','CellLabelColor','none');

%%

for gene=1:length(l)
    if ~isempty(find(strcmp(l{gene},chemokines)))
        ligands(gene,1)=find(strcmp(l{gene},chemokines));
        for clst = 1:12

             chemokines_map_Young(gene,1,clst) = mean(young{clst}(chemokines_idx(ligands(gene,1)),:),2); % For expression level analysis
             chemokines_map_Old(gene,1,clst) = mean(old{clst}(chemokines_idx(ligands(gene,1)),:),2); % For expression level analysis
             chemokines_map_h(gene,1,clst) = PValues{clst}(chemokines_idx(ligands(gene,1))); % For expression level analysis
            
        end
    else
        ligands(gene,1)=0;
        chemokines_map_Young(gene,1,1:12) = NaN;
        chemokines_map_Old(gene,1,1:12) = NaN;
        chemokines_map_h(gene,1,1:12) = NaN;
    end
    if ~isempty(find(strcmp(r{gene},chemokines)))
        receptors(gene,1)=find(strcmp(r{gene},chemokines));
        for clst = 1:12

             chemokines_map_Young(gene,2,clst) = mean(young{clst}(chemokines_idx(receptors(gene,1)),:),2); % For exoression level analysis   %Macrophages_age_changes{clst,1}
             chemokines_map_Old(gene,2,clst) = mean(old{clst}(chemokines_idx(receptors(gene,1)),:),2); % For exoression level analysis   %Macrophages_age_changes{clst,1}
             chemokines_map_h(gene,2,clst) = PValues{clst}(chemokines_idx(receptors(gene,1))); % For expression level analysis

        end
    else
        receptors(gene,1)=0;
        chemokines_map_Young(gene,2,1:12) = NaN;
        chemokines_map_Old(gene,2,1:12) = NaN;
        chemokines_map_h(gene,2,1:12) = NaN;
    end
    
end
%% FC
clear Chemokines_FC
clear Chemokines_score
clear Chemokines_sig
clear chemokines_final
clear x

Chemokines_FC = chemokines_map_Old./chemokines_map_Young;

for k = 1:size(Chemokines_FC,1) % Edge
    for clst=1:12 % Celltype i
        for j=1:12 % Celltype j
            if Chemokines_FC(k,1,clst)>= 2 && Chemokines_FC(k,2,j)>= 2 && Chemokines_FC(k,1,clst)< Inf && Chemokines_FC(k,2,j)< Inf
                Chemokines_score(clst,j,k) = 1;
            elseif Chemokines_FC(k,1,clst)<= .5 && Chemokines_FC(k,2,j)<= .5 && Chemokines_FC(k,1,clst)> 0 && Chemokines_FC(k,2,j)> 0
                Chemokines_score(clst,j,k) = -1;
            else
                Chemokines_score(clst,j,k) = 0;
            end
            
            if chemokines_map_h(k,1,clst)<= .05 && chemokines_map_h(k,2,j)<= .05
                Chemokines_sig(clst,j,k) = 1;
            else
                Chemokines_sig(clst,j,k) = 0;
            end
        end
    end
end

chemokines_final = Chemokines_sig.*Chemokines_score;

x=chemokines_final;
thresh = -1;
x(x > thresh) = 0;
x(x <=  thresh) = 1;

%% Exporting data to R
clear adjacency_list
l=0;
list_type = 1;
letters = {'A','B','C','D','E','F','G','H','I','J','K','L'};
for k=1:size(x,3)
    if ~isempty(find(x(:,:,k)))
        [row, col] = find(x(:,:,k));
        for num = 1:length(row)
            if list_type ==1
            adjacency_list{l+num,1} = [num2str(row(num)) ' - ' c{k+1,1}];
            adjacency_list{l+num,2} = [num2str(col(num)) ' - ' c{k+1,7}];
            adjacency_list{l+num,3} = 1;
            elseif list_type ==2
            adjacency_list{l+num,1} = c{k+1,1};
            adjacency_list{l+num,2} = c{k+1,7};
            adjacency_list{l+num,3} = 1;
            elseif list_type ==3
            adjacency_list{l+num,1} = cell_type_2clust{row(num)};
            adjacency_list{l+num,2} = cell_type_2clust{col(num)};
            adjacency_list{l+num,3} = 1;
            elseif list_type ==4
            adjacency_list{l+num,1} = [num2str(row(num)) ' - ' num2str(k)];
            adjacency_list{l+num,2} = [num2str(col(num)) ' - ' num2str(k)];
            adjacency_list{l+num,3} = 1;
            end
        end
        
        if k==1
            l=l+length(row)+1;
        else
            l=l+length(row);
        end
    end
end

%%

%%%%%%%%%%%%%%%%%%
%    Figure 5A   %
%%%%%%%%%%%%%%%%%%

clear l
clear r
clear c
clear cytokines_map
clear no_zero_clst
clear no_zero_cytokines
clear lig_flag
% Cytokines
[a,b,c] =xlsread('Cytokine_network_template.xlsx');

l = c(2:end,1);
r = c(2:end,7);

for clst=1:12
    for gene = 1:length(cytokines)
        FC = mean(old{clst}(cytokines_idx(gene),:),2)/mean(young{clst}(cytokines_idx(gene),:),2);
        if  PValues{clst}(cytokines_idx(gene))<.05 && FC>=2 && FC<Inf
            if -log10(PValues{clst}(cytokines_idx(gene)))<=3
                cytokines_map(clst,gene) = 1;
                logged(clst,gene) = -log10(PValues{clst}(cytokines_idx(gene)));
            elseif -log10(PValues{clst}(cytokines_idx(gene)))>3 && -log10(PValues{clst}(cytokines_idx(gene)))<=6
                cytokines_map(clst,gene) = 2;
                logged(clst,gene) = -log10(PValues{clst}(cytokines_idx(gene)));
            elseif -log10(PValues{clst}(cytokines_idx(gene)))>6
                cytokines_map(clst,gene) = 3;
                logged(clst,gene) = -log10(PValues{clst}(cytokines_idx(gene)));
                
            end
        elseif PValues{clst}(cytokines_idx(gene))<.05 && FC<=.5 && FC>0
            if -log10(PValues{clst}(cytokines_idx(gene)))<=3
                cytokines_map(clst,gene) = -1;
                logged(clst,gene) = log10(PValues{clst}(cytokines_idx(gene)));
                
            elseif -log10(PValues{clst}(cytokines_idx(gene)))>3 && -log10(PValues{clst}(cytokines_idx(gene)))<=6
                cytokines_map(clst,gene) = -2;
                logged(clst,gene) = log10(PValues{clst}(cytokines_idx(gene)));
                
            elseif -log10(PValues{clst}(cytokines_idx(gene)))>6 && -log10(PValues{clst}(cytokines_idx(gene)))<=9
                cytokines_map(clst,gene) = -3;
                logged(clst,gene) = log10(PValues{clst}(cytokines_idx(gene)));
                
            elseif -log10(PValues{clst}(cytokines_idx(gene)))>9
                cytokines_map(clst,gene) = -4;
                logged(clst,gene) = -log10(PValues{clst}(cytokines_idx(gene)));
                
            end
        else
            cytokines_map(clst,gene) = 0;
        end
    end
end

no_zero_clst_cyto = find(sum(cytokines_map~=0));
no_zero_cytokines = cytokines(no_zero_clst_cyto);

lig = 1;
rec = 1;
for clst=1:length(no_zero_cytokines)
   if sum(strcmp(no_zero_cytokines{clst},l))>0
    lig_flag_cyto(lig) = clst;
    lig=lig+1;
   elseif sum(strcmp(no_zero_cytokines{clst},r))>0
    rec_flag_cyto(rec) = clst;
    rec=rec+1;
   end
end
no_zero_cytokines_ordered = no_zero_cytokines([lig_flag_cyto rec_flag_cyto]);
no_zero_clst_ordered_cyto = no_zero_clst_cyto([lig_flag_cyto rec_flag_cyto]);

clr_map =[1 .2 .2; 1 .4 .4; 1 .6 .6; 1 .8 .8; 1 1 1; .9 1 .8; .8 1 .6; .7 1 .4];%; .6 1 .2;];
figure(26);
h1 = heatmap(cytokines_map(:,no_zero_clst_ordered_cyto),'Colormap',clr_map,'ColorbarVisible','off','CellLabelColor','none');
cdl = h1.XDisplayLabels; 
h1.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));
cdl = h1.YDisplayLabels; 
h1.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));
%%
clr_map = [1 .2 .2; 1 .4 .4; 1 .6 .6; 1 .8 .8; 1 1 1; .8 1 .6];
figure(30);
heatmap([chemokines_map(:,no_zero_clst_chmo(lig_flag_chmo)) cytokines_map(:,no_zero_clst_cyto(lig_flag_cyto))],'Colormap',clr_map,'ColorbarVisible','off','CellLabelColor','none');

clr_map =[1 .2 .2; 1 .4 .4; 1 .6 .6; 1 .8 .8; 1 1 1; .9 1 .8; .8 1 .6; .7 1 .4];%; .6 1 .2;];
figure(31);
heatmap([chemokines_map(:,no_zero_clst_chmo([16 rec_flag_chmo])) cytokines_map(:,no_zero_clst_cyto(rec_flag_cyto))],'Colormap',clr_map,'ColorbarVisible','off','CellLabelColor','none');

%%
clear ligands
for gene=1:length(l)
    if ~isempty(find(strcmp(l{gene},cytokines)))
        ligands(gene,1)=find(strcmp(l{gene},cytokines));

        for clst = 1:12

             cytokines_map_Young(gene,1,clst) = mean(young{clst}(cytokines_idx(ligands(gene,1)),:),2); % For expression level analysis
             cytokines_map_Old(gene,1,clst) = mean(old{clst}(cytokines_idx(ligands(gene,1)),:),2); % For expression level analysis
             cytokines_map_h(gene,1,clst) = PValues{clst}(cytokines_idx(ligands(gene,1))); % For expression level analysis

            
        end
    else
        ligands(gene,1)=0;
        cytokines_map_Young(gene,1,1:12) = NaN;
        cytokines_map_Old(gene,1,1:12) = NaN;
        cytokines_map_h(gene,1,1:12) = 1;
    end
    if ~isempty(find(strcmp(r{gene},cytokines)))
        receptors(gene,1)=find(strcmp(r{gene},cytokines));
        for clst = 1:12

             cytokines_map_Young(gene,2,clst) = mean(young{clst}(cytokines_idx(receptors(gene,1)),:),2); % For exoression level analysis   %Macrophages_age_changes{clst,1}
             cytokines_map_Old(gene,2,clst) = mean(old{clst}(cytokines_idx(receptors(gene,1)),:),2); % For exoression level analysis   %Macrophages_age_changes{clst,1}
             cytokines_map_h(gene,2,clst) = PValues{clst}(cytokines_idx(receptors(gene,1))); % For expression level analysis

        end
    else
        receptors(gene,1)=0;
        cytokines_map_Young(gene,2,1:12) = NaN;
        cytokines_map_Old(gene,2,1:12) = NaN;
        cytokines_map_h(gene,2,1:12) = 1;
    end
    
end
%% FC
clear Cytokines_FC
clear Cytokines_score
clear Cytokines_sig
clear chemokines_final
clear x

Cytokines_FC = cytokines_map_Old./cytokines_map_Young;

for k = 1:size(Cytokines_FC,1) % Edge
    for clst=1:12 % Celltype i
        for j=1:12 % Celltype j
            if Cytokines_FC(k,1,clst)>= 2 && Cytokines_FC(k,2,j)>= 2 && Cytokines_FC(k,1,clst)< Inf && Cytokines_FC(k,2,j)< Inf
                Cytokines_score(clst,j,k) = 1;
            elseif Cytokines_FC(k,1,clst)<= .5 && Cytokines_FC(k,2,j)<= .5 && Cytokines_FC(k,1,clst)> 0 && Cytokines_FC(k,2,j)> 0
                Cytokines_score(clst,j,k) = -1;
            else
                Cytokines_score(clst,j,k) = 0;
            end
            
            if cytokines_map_h(k,1,clst)<= .05 && cytokines_map_h(k,2,j)<= .05
                Cytokines_sig(clst,j,k) = 1;
            else
                Cytokines_sig(clst,j,k) = 0;
            end
        end
    end
end

cytokines_final = Cytokines_sig.*Cytokines_score;

x=cytokines_final;
thresh = -1;
x(x > thresh) = 0;
x(x <=  thresh) = 1;

%% Exporting data to R
clear adjacency_list
l=0;
list_type = 1;
letters = {'A','B','C','D','E','F','G','H','I','J','K','L'};
for k=1:size(x,3)
    if ~isempty(find(x(:,:,k)))
        [row, col] = find(x(:,:,k));
        for num = 1:length(row)
            if list_type ==1
            adjacency_list{l+num,1} = [num2str(row(num)) ' - ' c{k+1,1}];
            adjacency_list{l+num,2} = [num2str(col(num)) ' - ' c{k+1,7}];
            adjacency_list{l+num,3} = 1;
            elseif list_type ==2
            adjacency_list{l+num,1} = c{k+1,1};
            adjacency_list{l+num,2} = c{k+1,7};
            adjacency_list{l+num,3} = 1;
            elseif list_type ==3
            adjacency_list{l+num,1} = cell_type_2clust{row(num)};
            adjacency_list{l+num,2} = cell_type_2clust{col(num)};
            adjacency_list{l+num,3} = 1;
            end
        end
        
        if k==1
            l=l+length(row)+1;
        else
            l=l+length(row);
        end
    end
end

%% Setting significant genes based on diff (old-young)
 
% %%% https://www.mathworks.com/matlabcentral/fileexchange/45966-compare-two-proportions-chi-square %%%
% % Chi-square proportion test for significance testing
tbl = 1;
col=4;
for clst = 1:12
    for gene = 1:length(name_joint)
        pop_size_1 = length(poptype_10X_R{clst,1});
        pop_size_2 = length(poptype_10X_R{clst,2});
        N = [pop_size_1 pop_size_2];
        
        pos_1 = sum(data_both_expressed{1}(gene,poptype_10X_R{clst,1})>0);
        pos_2 = sum(data_both_expressed{2}(gene,poptype_10X_R{clst,2})>0);
        X = [pos_1 pos_2];
        
        fraction_young{clst}(gene,1) = 100*(pos_1/pop_size_1);
        fraction_old{clst}(gene,1) = 100*(pos_2/pop_size_2);
        fraction_diff{clst}(gene,1) = fraction_old{clst}(gene,1)-fraction_young{clst}(gene,1);
        
        if sum(X)~=0
            [h(gene,clst),p(gene,clst),chi2stat(gene,clst),df] = prop_test(X ,N, 'false');
        else
            h(gene,clst)=0;
            p(gene,clst)=1;
            chi2stat(gene,clst)=NaN;
        end
    FC_fractions{clst}(gene) = (pos_2/pop_size_2)/(pos_1/pop_size_1);   
    end
    
    up_ind{clst} = find(fraction_diff{clst}>=0);
    down_ind{clst} = find(fraction_diff{clst}<=0);
    
    up_genes_ind{clst} = intersect(up_ind{clst},find(h(:,clst)));
    down_genes_ind{clst} = intersect(down_ind{clst},find(h(:,clst)));
    
    [up_sort{clst},up_sort_ind{clst}]= sort(fraction_diff{clst}(up_genes_ind{clst}),'descend');
    [down_sort{clst},down_sort_ind{clst}]= sort(fraction_diff{clst}(down_genes_ind{clst}),'ascend');
    
    up_genes_name{clst} = name_joint(up_genes_ind{clst}(up_sort_ind{clst}));
    down_genes_name{clst} = name_joint(down_genes_ind{clst}(down_sort_ind{clst}));
    
    xplot = fraction_diff{clst};
    yplot = -log10(p(:,clst));
    Sigplot = [up_genes_ind{clst}; down_genes_ind{clst}];
    
    fraction_mat{clst,1} = up_genes_name{clst};
    fraction_mat{clst,2} = p(up_genes_ind{clst},clst);
    fraction_mat{clst,3} = fraction_diff{clst}(up_genes_ind{clst});
    fraction_mat{clst,4} = down_genes_name{clst};
    fraction_mat{clst,5} = p(down_genes_ind{clst},clst);
    fraction_mat{clst,6} = fraction_diff{clst}(down_genes_ind{clst});
end    

%%
SASP_genes = {'Ighm','Ccr2','Csf2ra','Clec4a2','Clec4a3','Ifngr1','Cxcr6','Csf1r'};
old_macs = poptype_10X_R{2,2};
young_macs = poptype_10X_R{2,1};

flag = 0;
for gene = 1:length(SASP_genes)
    ind_SASPgene(gene) = find(strcmp(SASP_genes{gene},name_joint));
    if ismember(ind_SASPgene(gene),up_genes_ind{2})
        flag=flag+1;
    end
end
flag

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Spplementary Figure 5     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f,dF] = ecdf(FrcDiff{2});

figure(27);
plot(dF,f,'LineWidth',2); hold on

for i=1:length(ind_SASPgene)
    [val(i),ind(i)] = min(abs(dF-FrcDiff{2}(ind_SASPgene(i))));
    pvl(i) = 1-f(ind(i));
    
end

scatter(dF(ind),f(ind),50,'rd','filled')
% text(dF(ind),f(ind)+0.05,SASP_genes,'FontSize',13)
ylim([0 1.05])
