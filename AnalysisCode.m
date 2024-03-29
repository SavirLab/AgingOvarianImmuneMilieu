%% Paper Graphs Script
close all
clc
clear

%% Loading Raw Data
DataLoading

%% Loading Suerat (R) parameters

[~,~,raw_cluster] = xlsread('Cluster_idnts.csv'); % Cluster ID for each cell
[~,~,raw_tsne] = xlsread('tSNE_coordinates.csv'); % tSNE coordinates
[~,~,raw_age] = xlsread('Age_vec.csv'); % Sample (age) ID for each cell

all_tSNE = cell2mat(raw_tsne);
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
    
    % We are using the following ss: no less than 200 and no more than 2500 features & less than 10% mitochondrial genes
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

%%

%%%%%%%%%%%%%%%%%%
%   Figure 1B    %
%%%%%%%%%%%%%%%%%%
num_clst = max(all_cluster);

figure(3);
for clst = 1:num_clst
scatter(all_tSNE(all_cluster==clst,1),all_tSNE(all_cluster==clst,2),5,color(clst,:),'o','filled'); hold on  %
end


%% Violin plot for cluster identification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Spplementary Figure 1     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    set(gca,'xtick',[1:14],'xticklabels',[])
    set(gca,'ytick',[0 top_value(cnt_gene)],'yticklabels',[]);
    xlim([0 14.5])
    ylabel(marker_name{cnt_gene})
end

%% SingleR classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Spplementary Figure 2     %
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

title('SingleR FINE classifiction results on total 10X data')

%% CD8- CD4- in dofferent tissues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Spplementary Figure 4     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%   Figure 1C    %
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
%   Figure 1C    %
%%%%%%%%%%%%%%%%%%

figure(7);
scatter(all_tSNE(age_num==2,1),all_tSNE(age_num==2,2),10,[1 .6 .8],'o','filled','MarkerFaceAlpha',.5); hold on %[1 .6 .2]
scatter(all_tSNE(age_num==1,1),all_tSNE(age_num==1,2),10,[.4 .7 1],'o','filled','MarkerFaceAlpha',.5); 
title('tSNE on joint non-zero genes - log(Data+1)+PCA15')
legend({'Young','Old'});

%% Calculating Fractions + Bar plot
%%%%%%%%%%%%%%%%%%
%   Figure 1D    %
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

        for clst=1:size(poptype_10X_R,1)
        clust_numbers(clst,age) = length(poptype_10X_R{clst,age});
        if isempty(fractions(clst,age))
            clust_numbers(clst,age) =0;
        end
    end
    clust_numbers(clst+1,age) = sum(age_num==age)-sum(clust_numbers(1:clst,age));
    
    
end
cell_type_2clust{13} = 'Others';

age_group = {'Young','Old'};
figure(8); hold on
b=bar(fractions(:,:));
b(1).FaceColor = [.4 .7 1];
b(2).FaceColor = [1 .6 .8];
set(gca,'xtick',1:13,'xticklabels',cell_type_2clust)
set(gca,'XTickLabelRotation',-30)
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels1 = string(poptype_num(:,1));
labels2 = string(poptype_num(:,2));

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Figure 1E or supplementary 5    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Violin Plot - Facs (Macs+all CD3+ cells):

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
cats = ({'Macrophages','CD3 lymphocytes'});

for age = 1:3
    for stage = 1:4
        Macs_mean(stage,age) = nanmean(cells_all(cells_all(:,5)==stage & cells_all(:,4)==age,1));
        lymph_mean(stage,age) = nanmean(cells_all(cells_all(:,5)==stage & cells_all(:,4)==age,2));
    end
end

figure(9);
violins_macs = violinplot(Macs_violin, cats,'ViolinColor',[.3 .6 0],'Width',0.1,'ShowMean',true); hold on 
violins_lymph = violinplot(Lymph_violin, cats, 'ViolinColor',[1 .7 .4],'Width',0.1,'ShowMean',true); 
% run to add cycle stages:
% plot(1:3,Macs_mean,'color',[.3 .6 0],'LineWidth',1.5)
% plot(1:3,lymph_mean,'color',[1 .7 .4],'LineWidth',1.5)
%

Young_sample = find(cells_all(:,4)==1);
Adult_sample = find(cells_all(:,4)==2);
Old_sample = find(cells_all(:,4)==3);

% Statistical analysis
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
%   Figure 1F    %
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

%% Finding DEGs

name_joint = gene_names_young(expressed_genes_joint); % Names of all expressed genes

% Calculating pVals including non-expressing cells
for clst = 1:size(poptype_10X_R,1)
    young{clst} = data_both_expressed{1}(:,poptype_10X_R{clst,1});
    old{clst} = data_both_expressed{2}(:,poptype_10X_R{clst,2});
    
    % significant DEGs (t-test)
    for gene = 1:length(expressed_genes_joint)
        [PValues{clst}(gene), TScores{clst}(gene)] = mattest(old{clst}(gene,:),young{clst}(gene,:),'VarType','unequal');
                
    end
    
    
    % Calculating Qvalue (FDR)
    [fdr,Q] = mafdr(PValues{clst});
    ind_01 = find(Q>0.1);
    Q(ind_01) = floor(100*Q(ind_01))/100;
    Qvalues(:, clst) = Q.';
    
    
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
    sig_up{clst} = intersect(intersect(find(PValues{clst}<=.05),find(Qvalues(:,clst)<=.1)),find(log2(expression_FC{clst})>=1));
    sig_up_name{clst} = name_joint(sig_up{clst});
    sig_up_FC{clst} = expression_FC{clst}(sig_up{clst});
    sig_up_pval{clst}(:,1) = PValues{clst}(sig_up{clst});
    sig_up_Qval{clst}(:,1) = Qvalues(sig_up{clst},clst);
    sig_down{clst} = intersect(intersect(find(PValues{clst}<=.05),find(Qvalues(:,clst)<=.1)),find(log2(expression_FC{clst})<=-1));
    sig_down_name{clst} = name_joint(sig_down{clst});
    sig_down_FC{clst} = expression_FC{clst}(sig_down{clst});
    sig_down_pval{clst}(:,1) = PValues{clst}(sig_down{clst});
    sig_down_Qval{clst}(:,1) = Qvalues(sig_down{clst},clst);
    
    % expression proportion change (chi-square)
    expression_fraction{clst,1} = 100*(sum(young{clst}>0,2)./size(young{clst},2));
    expression_fraction{clst,2} = 100*(sum(old{clst}>0,2)./size(old{clst},2));
    
    clst
end

%%
%%%%%%%%%%%%%%%%%%
%    Figure 3A    %
%%%%%%%%%%%%%%%%%%

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

red_dots{1} = [ind_Ifngr1 ind_Tgfbr1 ind_Cxcl2  ind_Il1a...
    ind_Il1b ind_Csf2rb];
red_dots_names{1} = {'Ifngr1','Tgfbr1','Cxcl2','Il1a','Il1b',...
    'Csf2rb'};

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
ind_Ccl4 = find(strcmp('Ccl4',name_joint)); % Macs
ind_Inhba = find(strcmp('Inhba',name_joint)); % Macs

red_dots{2} = [ind_Ifngr1 ind_Cxcl2 ind_Cxcl10 ind_Il2ra ...
           ind_Tnf ind_Ccl4 ind_Inhba];
red_dots_names{2} = {'Ifngr1','Cxcl2','Cxcl10','Il2ra','Tnf',...
    'Ccl4','Inhba'};

ind_Cd81 = find(strcmp('Cd81',name_joint)); % DC
ind_F2r = find(strcmp('F2r',name_joint)); % DC
ind_Ccr5 = find(strcmp('Ccr5',name_joint)); % DC, ILC1
ind_Cxcl16 = find(strcmp('Cxcl16',name_joint)); % DC
ind_Stat1 = find(strcmp('Stat1',name_joint)); % DC,NK,ILC1
ind_Ccl3 = find(strcmp('Ccl3',name_joint)); % DC,NK,ILC1

red_dots{3} = [ind_Cxcl2 ind_Ccr5 ind_Cxcl16 ind_Ccl4 ind_Ccl3];
red_dots_names{3} = {'Cxcl2','Ccr5','Cxcl16','Ccl4','Ccl3'};

ind_Stat4 = find(strcmp('Stat4',name_joint)); % NK
ind_Icam1 = find(strcmp('Icam1',name_joint)); % NK
ind_Ccl3 = find(strcmp('Ccl3',name_joint)); % NK

red_dots{4} = [ind_Cxcl2 ind_Ccl3];
red_dots_names{4} = {'Cxcl2','Ccl3'};

ind_Cxcr4 = find(strcmp('Cxcr4',name_joint)); % NK
ind_Ccl7 = find(strcmp('Ccl7',name_joint)); % NK
ind_Ccl9 = find(strcmp('Ccl9',name_joint)); % NK

red_dots{5} = [ind_Cxcl2 ind_Il1b ind_Ccr5 ind_Tgfbr1 ind_Cxcr4 ind_Ccl7];
red_dots_names{5} = {'Cxcl2','Il1b','Ccr5','Tgfbr1','Cxcr4','Ccl7'};

ind_Ccl5 = find(strcmp('Ccl5',name_joint)); % NKT, ILC2
ind_Fcer1g = find(strcmp('Fcer1g',name_joint)); % NKT
ind_Mmp9 = find(strcmp('Mmp9',name_joint)); % NKT
ind_Il2 = find(strcmp('Il2',name_joint)); % NKT

red_dots{6} = [ind_Cxcl2 ind_Il1b ind_Ccl5];
red_dots_names{6} = {'Cxcl2','Il1b','Ccl5'};

ind_Ly6c2 = find(strcmp('Ly6c2',name_joint)); % CD8, CD4
ind_Ccr2 = find(strcmp('Ccr2',name_joint)); % NKT, ILC2

red_dots{7} = [ind_Cxcl2 ind_Il1b ];
red_dots_names{7} = {'Cxcl2','Il1b'};

ind_Icam2 = find(strcmp('Icam2',name_joint)); % CD4, ILC2
ind_Ccr6 = find(strcmp('Ccr6',name_joint)); % NKT, ILC2

red_dots{8} = [ind_Cxcl2];
red_dots_names{8} = {'Cxcl2'};

ind_Mmp16 = find(strcmp('Mmp16',name_joint)); % DNT
ind_Ccl4 = find(strcmp('Ccl4',name_joint)); % DNT,B
ind_Cxcl10 = find(strcmp('Cxcl10',name_joint)); % DNT,B
ind_Il1r1 = find(strcmp('Il1r1',name_joint)); % DNT,B
ind_Tgfbr2 = find(strcmp('Tgfbr2',name_joint)); % DNT,B

red_dots{9} = [ind_Cxcl2 ind_Ccl4 ind_Ccl3 ind_Cxcl10 ind_Il1b...
    ind_Ccr6 ind_Il1r1 ind_Tgfbr2];
red_dots_names{9} = {'Cxcl2','Ccl4','Ccl3','Cxcl10','Il1b',...
    'Ccr6','Il1r1','Tgfbr2'};

ind_Iglc2 = find(strcmp('Iglc2',name_joint)); % B

red_dots{10} = [ind_Cxcl2];
red_dots_names{10} = {'Cxcl2'};

ind_Il10rb = find(strcmp('Il10rb',name_joint)); % ILC2
ind_Ccr2 = find(strcmp('Ccr2',name_joint)); % ILC2

red_dots{11} = [ind_Ccr2 ind_Ccl5];
red_dots_names{11} = {'Ccr2','Ccl5'};

ind_Xcl1 = find(strcmp('Xcl1',name_joint)); % ILC3
ind_Il7r = find(strcmp('Il7r',name_joint)); % ILC3
ind_Ccr6 = find(strcmp('Ccr6',name_joint)); % ILC3
ind_Prdx2 = find(strcmp('Prdx2',name_joint)); % ILC3

red_dots{12} = [ind_Il7r ind_Ccr6];
red_dots_names{12} = {'Il7r','Ccr6'};

%%
close all
q_line_01 = [0.005 0.022 0.066 0.005 0.040 0.013 0.018 0.006 0.057 0.0007 0.0075 0.0008];
for clst = 1:size(poptype_10X_R,1)
   x = log2(expression_FC{clst});
   y = -log10(PValues{clst});
   z=-log10(q_line_01(clst));
   
   figure(11+clst);
   scatter(x,y,8,[.6 .6 .6],'o','filled'); hold on
   scatter(x(red_dots{clst}),y(red_dots{clst})',10,'ro','filled');
   text(x(red_dots{clst})+0.1,y(red_dots{clst})+0.1,red_dots_names{clst},'FontSize',13)
   line([1 1],[min(y)-0.5 max(y)+0.5],'LineStyle','--')
   line([-1 -1],[min(y)-0.5 max(y)+0.5],'LineStyle','--')
   line([min(x(x>-Inf))-0.5 max(x(x<Inf))+0.5],[0.585 0.585],'LineStyle','--')
   line([min(x(x>-Inf))-0.5 max(x(x<Inf))+0.5],[z z],'color','red','LineStyle','--')
   ylim([0 max(y)+0.5])
   xlim([-6 6])
   title(cell_type_2clust{clst})
end
    
%% 
%%%%%%%%%%%%%%%%%%
%    Figure 3B   %
%%%%%%%%%%%%%%%%%%

[~,txt_REVIGO_common_down,raw_REVIGO_common_down] = xlsread([pwd '\REVIGO_Down_Common.xlsx']);

close all
figure(24)
for i=2:size(raw_REVIGO_common_down,1)
    if isnumeric(raw_REVIGO_common_down{i,8})
        scatter(raw_REVIGO_common_down{i,8},raw_REVIGO_common_down{i,9},2.6^raw_REVIGO_common_down{i,4},raw_REVIGO_common_down{i,3},...
            'o','filled','MarkerEdgeColor',[0 0 0]); hold on
    end
end
h = gca;
h.XAxis.Visible = 'off'; h.YAxis.Visible = 'off';


%% Finding Cytokines and Chemokines 
chemokines_prefix = {'Ccl','Ccrl','Cxcl','Cx3cl','Xcl','Xcr','Ackr','Ccr','Cxcr','Cx3cr'};
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
close all
%%
%%%%%%%%%%%%%%%%%%
%    Figure 4A   %
%%%%%%%%%%%%%%%%%%
qnum = Qvalues;

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
clear chemokines_map

for clst=1:12
    for gene = 1:length(chemokines)
        FC = mean(old{clst}(chemokines_idx(gene),:),2)/mean(young{clst}(chemokines_idx(gene),:),2);
       if  PValues{clst}(chemokines_idx(gene))<.05 && FC>=2 && FC<Inf
           if -log10(qnum(chemokines_idx(gene),clst))>=1 && -log10(qnum(chemokines_idx(gene),clst))<2 %-log10(PValues{clst}(chemokines_idx(gene)))<=3
               chemokines_map(clst,gene) = 1;
           elseif -log10(qnum(chemokines_idx(gene),clst))>=2 && -log10(qnum(chemokines_idx(gene),clst))<3 %-log10(PValues{clst}(chemokines_idx(gene)))>3 && -log10(PValues{clst}(chemokines_idx(gene)))<=6
               chemokines_map(clst,gene) = 2;
           elseif -log10(qnum(chemokines_idx(gene),clst))>=3 %-log10(PValues{clst}(chemokines_idx(gene)))>6
               chemokines_map(clst,gene) = 3;
           end
       elseif PValues{clst}(chemokines_idx(gene))<.05 && FC<=.5 && FC>0
           if -log10(qnum(chemokines_idx(gene),clst))>=1 && -log10(qnum(chemokines_idx(gene),clst))<2 %-log10(PValues{clst}(chemokines_idx(gene)))<=3
               chemokines_map(clst,gene) = -1;
           elseif -log10(qnum(chemokines_idx(gene),clst))>=2 && -log10(qnum(chemokines_idx(gene),clst))<3 %-log10(PValues{clst}(chemokines_idx(gene)))>3 && -log10(PValues{clst}(chemokines_idx(gene)))<=6
               chemokines_map(clst,gene) = -2;
           elseif -log10(qnum(chemokines_idx(gene),clst))>=3 && -log10(qnum(chemokines_idx(gene),clst))<4 %-log10(PValues{clst}(chemokines_idx(gene)))>6 && -log10(PValues{clst}(chemokines_idx(gene)))<=9
               chemokines_map(clst,gene) = -3;
           elseif -log10(qnum(chemokines_idx(gene),clst))>=4 %-log10(PValues{clst}(chemokines_idx(gene)))>9
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
lig_chemo = lig-1;
rec_chemo = rec-1;

no_zero_chemokines_ordered = no_zero_chemokines([lig_flag_chmo(1:lig_chemo) 16 rec_flag_chmo]);
no_zero_clst_ordered_chmo = no_zero_clst_chmo([lig_flag_chmo(1:lig_chemo) 16 rec_flag_chmo]);


clr_map = [1 .2 .2; 1 .4 .4; 1 .6 .6; 1 .8 .8; 1 1 1; .8 1 .6]; 
% figure(25);
% heatmap(chemokines_map(:,no_zero_clst_ordered_chmo),'Colormap',clr_map,'ColorbarVisible','off','CellLabelColor','none');

%% Creating table for R - chord graph, chemokines

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
%    Figure 4A   %
%%%%%%%%%%%%%%%%%%

%Cytokines
clear l
clear r
clear c
clear cytokines_map
clear no_zero_clst
clear no_zero_cytokines
clear lig_flag

[a,b,c] =xlsread('Cytokine_network_template.xlsx');

l = c(2:end,1);
r = c(2:end,7);

for clst=1:12
    for gene = 1:length(cytokines)
        FC = mean(old{clst}(cytokines_idx(gene),:),2)/mean(young{clst}(cytokines_idx(gene),:),2);
        if  PValues{clst}(cytokines_idx(gene))<.05 && FC>=2 && FC<Inf
            if -log10(qnum(cytokines_idx(gene),clst))>=1 && -log10(qnum(cytokines_idx(gene),clst))<2 %-log10(PValues{clst}(cytokines_idx(gene)))<=3
                cytokines_map(clst,gene) = 1;
                logged(clst,gene) = -log10(PValues{clst}(cytokines_idx(gene)));
            elseif -log10(qnum(cytokines_idx(gene),clst))>=2 && -log10(qnum(cytokines_idx(gene),clst))<3 %-log10(PValues{clst}(cytokines_idx(gene)))>3 && -log10(PValues{clst}(cytokines_idx(gene)))<=6
                cytokines_map(clst,gene) = 2;
                logged(clst,gene) = -log10(PValues{clst}(cytokines_idx(gene)));
            elseif -log10(qnum(cytokines_idx(gene),clst))>=3 %-log10(PValues{clst}(cytokines_idx(gene)))>6
                cytokines_map(clst,gene) = 3;
                logged(clst,gene) = -log10(PValues{clst}(cytokines_idx(gene)));
                
            end
        elseif PValues{clst}(cytokines_idx(gene))<.05 && FC<=.5 && FC>0
            if -log10(qnum(cytokines_idx(gene),clst))>=1 && -log10(qnum(cytokines_idx(gene),clst))<2 %-log10(PValues{clst}(cytokines_idx(gene)))<=3
                cytokines_map(clst,gene) = -1;
                logged(clst,gene) = log10(PValues{clst}(cytokines_idx(gene)));
                
            elseif -log10(qnum(cytokines_idx(gene),clst))>=2 && -log10(qnum(cytokines_idx(gene),clst))<3 %-log10(PValues{clst}(cytokines_idx(gene)))>3 && -log10(PValues{clst}(cytokines_idx(gene)))<=6
                cytokines_map(clst,gene) = -2;
                logged(clst,gene) = log10(PValues{clst}(cytokines_idx(gene)));
                
            elseif -log10(qnum(cytokines_idx(gene),clst))>=3 && -log10(qnum(cytokines_idx(gene),clst))<4 %-log10(PValues{clst}(cytokines_idx(gene)))>6 && -log10(PValues{clst}(cytokines_idx(gene)))<=9
                cytokines_map(clst,gene) = -3;
                logged(clst,gene) = log10(PValues{clst}(cytokines_idx(gene)));
                
            elseif -log10(qnum(cytokines_idx(gene),clst))>=4 %-log10(PValues{clst}(cytokines_idx(gene)))>9
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
lig_cyto = lig-1;

no_zero_cytokines_ordered = no_zero_cytokines([lig_flag_cyto rec_flag_cyto]);
no_zero_clst_ordered_cyto = no_zero_clst_cyto([lig_flag_cyto rec_flag_cyto]);

clr_map =[1 .2 .2; 1 .4 .4; 1 .6 .6; 1 .8 .8; 1 1 1; .9 1 .8; .8 1 .6; .7 1 .4];
% figure(26);
% h1 = heatmap(cytokines_map(:,no_zero_clst_ordered_cyto),'Colormap',clr_map,'ColorbarVisible','off','CellLabelColor','none');
% cdl = h1.XDisplayLabels; 
% h1.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));
% cdl = h1.YDisplayLabels; 
% h1.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));

%% Figure 4A
cyto_map = cytokines_map(:,no_zero_clst_ordered_cyto);
chmo_map = chemokines_map(:,no_zero_clst_ordered_chmo);

clr_map =[1 .2 .2; 1 .4 .4; 1 .6 .6; 1 .8 .8; 1 1 1; .9 1 .8; .8 1 .6; .7 1 .4];

figure(30);
heatmap([chmo_map(:,1:lig_chemo) cyto_map(:,1:lig_cyto)],'Colormap',clr_map,'ColorbarVisible','off','CellLabelColor','none');

figure(31);
heatmap([chmo_map(:,lig_chemo+1:end) cyto_map(:,lig_cyto+1:end)],'Colormap',clr_map,'ColorbarVisible','off','CellLabelColor','none');

%% Creating table for R - chord graph, cytokines
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
        
%         if sum(X)~=0
%             [h(gene,clst),p(gene,clst),chi2stat(gene,clst),df] = prop_test(X ,N, 'false');
%         else
%             h(gene,clst)=0;
%             p(gene,clst)=1;
%             chi2stat(gene,clst)=NaN;
%         end
    FC_fractions{clst}(gene) = (pos_2/pop_size_2)/(pos_1/pop_size_1);   
    end
    
%     up_ind{clst} = find(fraction_diff{clst}>=0);
%     down_ind{clst} = find(fraction_diff{clst}<=0);
    
%     up_genes_ind{clst} = intersect(up_ind{clst},find(h(:,clst)));
%     down_genes_ind{clst} = intersect(down_ind{clst},find(h(:,clst)));
%     
%     [up_sort{clst},up_sort_ind{clst}]= sort(fraction_diff{clst}(up_genes_ind{clst}),'descend');
%     [down_sort{clst},down_sort_ind{clst}]= sort(fraction_diff{clst}(down_genes_ind{clst}),'ascend');
%     
%     up_genes_name{clst} = name_joint(up_genes_ind{clst}(up_sort_ind{clst}));
%     down_genes_name{clst} = name_joint(down_genes_ind{clst}(down_sort_ind{clst}));
    
%     xplot = fraction_diff{clst};
%     yplot = -log10(p(:,clst));
%     Sigplot = [up_genes_ind{clst}; down_genes_ind{clst}];
    
%     fraction_mat{clst,1} = up_genes_name{clst};
%     fraction_mat{clst,2} = p(up_genes_ind{clst},clst);
%     fraction_mat{clst,3} = fraction_diff{clst}(up_genes_ind{clst});
%     fraction_mat{clst,4} = down_genes_name{clst};
%     fraction_mat{clst,5} = p(down_genes_ind{clst},clst);
%     fraction_mat{clst,6} = fraction_diff{clst}(down_genes_ind{clst});
end    


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SASPS Figure 5     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SASP_genes_combined = {'IL1R1'...
,'IL1R2'...
,'IL6RA'...
,'IL6ST'...
,'CXCR1'...
,'CXCR2'...
,'TGFBR1'...
,'TGFBR2'...
,'LRP1'...
,'TMEM219'...
,'Plaur'...
,'CCR5'...
,'CCR1'...
,'CCR3'...
,'CCR6'...
,'CD74'...
,'Ighm'...
,'Ccr2'...
,'Csf2ra'...
,'Clec4a2'...
,'Clec4a3'...
,'Ifngr1'...
,'Cxcr6'...
,'Csf1r'...
};

for gene = 1:length(SASP_genes_combined)
    if ~isempty(find(strcmpi(SASP_genes_combined{gene},name_joint)))
    ind_SASPgene_combined(gene) = find(strcmpi(SASP_genes_combined{gene},name_joint));
    else
        ind_SASPgene_combined(gene) = nan;
    end
end


close all
clear val ind
[f,dF] = ecdf(fraction_diff{2});

figure(27);
plot(dF,f,'LineWidth',2); hold on

for i=1:length(ind_SASPgene_combined)
    [val(i),ind(i)] = min(abs(dF-fraction_diff{2}(ind_SASPgene_combined(i))));
    pvl(i) = 1-f(ind(i));
    
end

% scatter(dF(ind),f(ind),50,'rd','filled')
% ylim([0 1.05])
%
% for i=1:length(ind_SASPgene_combined)
%     [val_2(i),ind_2(i)] = min(abs(dF-FrcDiff{2}(ind_SASPgene_combined(i))));
%     pvl_2(i) = 1-f(ind_2(i));
%
% end

scatter(dF(ind),f(ind),50,'gd','filled')
text(dF(ind),f(ind)+0.1*rand(24,1),SASP_genes_combined')
% 
% for k =1
%     rand_ind = randperm(18110);
%     rand_ind = rand_ind(1:24)
%     
%     for i=1:length(rand_ind)
%         [val(i),ind(i)] = min(abs(dF-FrcDiff{9}(rand_ind(i))));
%         pvl(i) = 1-f(ind(i));
%         
%         scatter(dF(ind),f(ind),50,'kd','filled')
%         
%     end
% end
%% random cdf 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SASPS Figure 5  - Supp Fig 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear df f h p n_crit
n_rep = 100;

data_mphs = log2(expression_FC{2});

[f,dF,flo,fup] = ecdf(data_mphs);


%data_mphs = log2(expression_FC{2});

    [h_crit,p_crit] = kstest2(data_mphs(ind_SASPgene_combined),data_mphs);

    [val, ind] = min(abs(f-0.99))

df_crit = dF(ind);


tic

for k = 1:n_rep
    
    clear pvl
    rand_ind = randperm(18110);
    rand_ind = rand_ind(1:24);
    
    [f,dF] = ecdf(data_mphs(rand_ind));

    figure(1)

    plot(dF,f,'LineWidth',0.5,'color',[0.8 0.8 0.8]); hold on
    
    [h(k),p(k)] = kstest2(data_mphs(rand_ind),data_mphs);
    
    n_crit(k) = length(find(data_mphs(rand_ind)>df_crit));

end

toc
for i = 0:max(n_crit)
    
    fdr(i+1) = sum(n_crit==i)/n_rep;
end

% figure
% bar([0:max(n_crit)],log10(fdr))



[f,dF] = ecdf(data_mphs(ind_SASPgene_combined));
    plot(dF,f,'LineWidth',0.5,'color',[0.5 0 0]); hold on
    
[f,dF,flo,fup] = ecdf(data_mphs);

    plot(dF,f,'LineWidth',0.5,'color',[0 0 0.5]); hold on

    
    
    
[val, ind] = min(abs(f-0.99))

df_crit = dF(ind);

ind_up = find(data_mphs>df_crit);
name_up = name_joint(ind_up);

length(find(p<=p_crit))/n_rep

 Set_fig_YS(figure(1),18,18,18);
 box off
 %% Check the effect of SAPS on all cell types
 
 
for cnt_cluster = 1:12
    
data_clust = fraction_diff{cnt_cluster};

[f,dF,flo,fup] = ecdf(data_clust);
    [h_crit(cnt_cluster),p_crit(cnt_cluster)] = kstest2(data_clust(ind_SASPgene_combined),data_clust);


%data_mphs = log2(expression_FC{2});

end


%%
%%%%%%%%%%%%%%%%
%%%compare batch correcytion slides.
%%%%%%%%%%%%%%%%%%

clear  age_index data_wo_bc data_w_bc
close all

[num,txt,raw] = xlsread('C:\Users\SavirLab\Technion\Yoni Savir - TalBenYakov\Paper submission\elife\review\GitHub files\age_factor.csv');
txt = txt(2:end)

age_index = zeros(size(txt));

age_index(find(strcmp(txt,'Old')))=1;

data_wo_bc = zeros([length(num),3]);
data_w_bc = zeros([length(num),3]);

% with bc
[num,txt,raw] = xlsread('C:\Users\SavirLab\Technion\Yoni Savir - TalBenYakov\Paper submission\elife\review\GitHub files\tsnecorr_w_batch.csv');
data_w_bc(:,1:2) = num;
[num,txt,raw] = xlsread('C:\Users\SavirLab\Technion\Yoni Savir - TalBenYakov\Paper submission\elife\review\GitHub files\tsnecorr_w_batch_clusters.csv');
data_w_bc(:,3) = num(:,2);

% without bc
[num,txt,raw] = xlsread('C:\Users\SavirLab\Technion\Yoni Savir - TalBenYakov\Paper submission\elife\review\GitHub files\tsnecorr_wo_batch.csv');
% data_wo_bc(:,1:2) = num; % take the "new" tsne coordinates from the R

data_wo_bc(:,1:2) = all_tSNE; % take Tal's original tsne coordinaties from the paper.

[num,txt,raw] = xlsread('C:\Users\SavirLab\Technion\Yoni Savir - TalBenYakov\Paper submission\elife\review\GitHub files\tsnecorr_wo_batch_clusters.csv');
data_wo_bc(:,3) = num(:,2);
%

num_clust = 16;



%

figure(1)

subplot(1,2,1);hold on
scatter(data_wo_bc(find(age_index==0),1),data_wo_bc(find(age_index==0),2),10,[.4 .7 1],'o','filled','MarkerFaceAlpha',.5); 
scatter(data_wo_bc(age_index==1,1),data_wo_bc(age_index==1,2),10,[1 .6 .8],'o','filled','MarkerFaceAlpha',.5); 
%    legend({'Young','Old'});
   title('w/o bc')
   
   subplot(1,2,2);hold on
scatter(data_w_bc(find(age_index==0),1),data_w_bc(find(age_index==0),2),10,[.4 .7 1],'o','filled','MarkerFaceAlpha',.5); 
scatter(data_w_bc(age_index==1,1),data_w_bc(age_index==1,2),10,[1 .6 .8],'o','filled','MarkerFaceAlpha',.5); 
%    legend({'Young','Old'});
      title('w bc')

for cnt_cluster = 1:num_clust
    
    ind_clust = find(data_wo_bc(:,3)==cnt_cluster-1);
    
    
    subplot(1,2,1)
   scatter(mean(data_wo_bc(ind_clust,1)),mean(data_wo_bc(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 

    text(mean(data_wo_bc(ind_clust,1)),mean(data_wo_bc(ind_clust,2)),num2str((cnt_cluster)))
    
    subplot(1,2,2)
    
       scatter(mean(data_w_bc(ind_clust,1)),mean(data_w_bc(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 

    text(mean(data_w_bc(ind_clust,1)),mean(data_w_bc(ind_clust,2)),num2str((cnt_cluster)))
    
    
    
end
%%

% remove 14 and 15
ind_clust = find(data_wo_bc(:,3)==13);
data_wo_bc(ind_clust,1:2) = nan;
data_w_bc(ind_clust,1:2) = nan;

ind_clust = find(data_wo_bc(:,3)==14);
data_wo_bc(ind_clust,1:2) = nan;
data_w_bc(ind_clust,1:2) = nan;


% 2 = 6
data_w_bc((data_w_bc(:,3)==2-1),3)=6-1;
% 11=9
data_w_bc((data_w_bc(:,3)==11-1),3)=9-1;

%16=1
data_w_bc((data_w_bc(:,3)==16-1),3)=1-1;

figure(2)

subplot(1,2,1);hold on
scatter(data_wo_bc(find(age_index==0),1),data_wo_bc(find(age_index==0),2),10,[.4 .7 1],'o','filled','MarkerFaceAlpha',.5); 
scatter(data_wo_bc(age_index==1,1),data_wo_bc(age_index==1,2),10,[1 .6 .8],'o','filled','MarkerFaceAlpha',.5); 
%    legend({'Young','Old'});
%    title('w/o bc')
   
xlim([-50 50]);ylim([-60 60]);
axis square
Set_fig_YS(figure(2),18,18,18)
set(gca,'visible','off')

   subplot(1,2,2);hold on
scatter(data_w_bc(find(age_index==0),1),data_w_bc(find(age_index==0),2),10,[.4 .7 1],'o','filled','MarkerFaceAlpha',.5); 
scatter(data_w_bc(age_index==1,1),data_w_bc(age_index==1,2),10,[1 .6 .8],'o','filled','MarkerFaceAlpha',.5); 
%    legend({'Young','Old'});
%       title('w bc')
xlim([-50 50]);ylim([-60 60]);
axis square
Set_fig_YS(figure(2),18,18,18)
set(gca,'visible','off')

% for cnt_cluster = 1:num_clust
%     
%     ind_clust = find(data_wo_bc(:,3)==cnt_cluster-1);
%     
%     
%     subplot(1,2,1)
%    scatter(nanmean(data_wo_bc(ind_clust,1)),nanmean(data_wo_bc(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 
% 
%     text(nanmean(data_wo_bc(ind_clust,1)),nanmean(data_wo_bc(ind_clust,2)),num2str((cnt_cluster-1)))
%     
%     
%         ind_clust = find(data_w_bc(:,3)==cnt_cluster-1);
% 
%     subplot(1,2,2)
%     
%        scatter(nanmean(data_w_bc(ind_clust,1)),nanmean(data_w_bc(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 
% 
%     text(nanmean(data_w_bc(ind_clust,1)),nanmean(data_w_bc(ind_clust,2)),num2str((cnt_cluster-1)))
%     
%     
%     
% end



%


%
figure(3)
for cnt_cluster = 1:num_clust
    
    ind_clust = find(all_cluster==cnt_cluster-1);
    
    
    subplot(1,2,1);hold on
       scatter((data_wo_bc(ind_clust,1)),(data_wo_bc(ind_clust,2)),10,'o','filled','MarkerFaceAlpha',.5); 

   scatter(nanmean(data_wo_bc(ind_clust,1)),nanmean(data_wo_bc(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 

    text(nanmean(data_wo_bc(ind_clust,1)),nanmean(data_wo_bc(ind_clust,2)),num2str((cnt_cluster-1)))
    
    
    
    ind_clust = find(data_w_bc(:,3)==cnt_cluster-1);

    
    subplot(1,2,2);hold on
    
           scatter((data_w_bc(ind_clust,1)),(data_w_bc(ind_clust,2)),10,'o','filled','MarkerFaceAlpha',.5); 

       scatter(nanmean(data_w_bc(ind_clust,1)),nanmean(data_w_bc(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 

    text(nanmean(data_w_bc(ind_clust,1)),nanmean(data_w_bc(ind_clust,2)),num2str((cnt_cluster-1)))
    
    
    
end



%% allign w and w/o cluster numbers
num_clust = 16;
data_w_bc_aligned = data_w_bc;

  %   0     2     3     4     5     6     7     8     9    11    12    13    14

align_cluster =  [1 4 5 6 2 7 9 14 11 13 12];

% check allign

cluster_vec = unique(data_wo_bc(:,3))';

for cnt_cluster = 1:length(cluster_vec)
    
    ind_cluster = find(data_w_bc(:,3)==cluster_vec(cnt_cluster));
    
    align_cluster_auto(cnt_cluster) = round(nanmedian(all_cluster(ind_cluster)));
    
    
    
if ~isempty(ind_cluster)
    data_w_bc_aligned(ind_cluster,3) = align_cluster_auto(cnt_cluster);
end
    
end
    
%%
close all

color_batch = [1 .2 .6; 1 .4 .4; 1 .6 .2; .75 .75 .75;...
    .6 .6 0; .3 .6 0; .6 .3 0; 0 .6 .3;...
    0 .6 .6; 0 .8 .8; 0 0 0; .4 .7 1;...
    .6 .6 1; 1 .6 1;];

figure(4)
for cnt_cluster = 1:num_clust
    
    ind_clust = find(all_cluster==cnt_cluster-1);
    
      
    subplot(1,2,1);hold on
       scatter((data_wo_bc(ind_clust,1)),(data_wo_bc(ind_clust,2)),10,color_batch(cnt_cluster,:),'o','filled','MarkerFaceAlpha',.5); 

%    scatter(nanmean(data_wo_bc(ind_clust,1)),nanmean(data_wo_bc(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 

%     text(nanmean(data_wo_bc(ind_clust,1)),nanmean(data_wo_bc(ind_clust,2)),num2str((cnt_cluster)))
    
    
    
    ind_clust = find(data_w_bc_aligned(:,3)==cnt_cluster-1);

    
    subplot(1,2,2);hold on
    
           scatter((data_w_bc_aligned(ind_clust,1)),(data_w_bc_aligned(ind_clust,2)),10,color_batch(cnt_cluster,:),'o','filled','MarkerFaceAlpha',.5); 

%        scatter(nanmean(data_w_bc(ind_clust,1)),nanmean(data_w_bc_aligned(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 

%     text(nanmean(data_w_bc_aligned(ind_clust,1)),nanmean(data_w_bc_aligned(ind_clust,2)),num2str((cnt_cluster)))
    
    
    
end
figure(4)
   subplot(1,2,1)
xlim([-50 50]);ylim([-60 60]);
axis square
Set_fig_YS(figure(4),18,18,18)
set(gca,'visible','off')

   subplot(1,2,2)
xlim([-50 50]);ylim([-60 60]);
axis square
Set_fig_YS(figure(4),18,18,18)
set(gca,'visible','off')

%%
close all

color = [1 .4 .4; 1 .6 .2; .75 .75 .75; .6 .6 0;...
    .3 .6 0; .6 .3 0; 0 .6 .3; 0 .6 .6;...
    0 .8 .8; .2 .6 1; .4 .7 1; .6 .6 1;...
    1 .6 1; 1 .2 .6];

color_batch = [1 .2 .6; 1 .4 .4; 1 .6 .2; .75 .75 .75;...
    .6 .6 0; .3 .6 0; .6 .3 0; 0 .6 .3;...
    0 .6 .6; 0 .8 .8; 0 0 0; .4 .7 1;...
    .6 .6 1; 1 .6 1;];

figure(5)
for cnt_cluster = 1:num_clust
    
    ind_clust = find(all_cluster==cnt_cluster);
    
    
    subplot(1,2,1);hold on
       scatter((data_wo_bc(ind_clust,1)),(data_wo_bc(ind_clust,2)),10,color(cnt_cluster,:),'o','filled','MarkerFaceAlpha',.5); 

%    scatter(nanmean(data_wo_bc(ind_clust,1)),nanmean(data_wo_bc(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 

%     text(nanmean(data_wo_bc(ind_clust,1)),nanmean(data_wo_bc(ind_clust,2)),num2str((cnt_cluster)))
    
    
    
%     ind_clust = find(data_w_bc_aligned(:,3)==cnt_cluster-1);

    
    subplot(1,2,2);hold on
    
           scatter((data_w_bc_aligned(ind_clust,1)),(data_w_bc_aligned(ind_clust,2)),10,color(cnt_cluster,:),'o','filled','MarkerFaceAlpha',.5); 

%        scatter(nanmean(data_w_bc(ind_clust,1)),nanmean(data_w_bc_aligned(ind_clust,2)),10,'k','o','filled','MarkerFaceAlpha',.5); 

    text(nanmean(data_w_bc_aligned(ind_clust,1)),nanmean(data_w_bc_aligned(ind_clust,2)),num2str((cnt_cluster)))
    
    
    
end
figure(5)
   subplot(1,2,1)
xlim([-50 50]);ylim([-60 60]);
axis square
Set_fig_YS(figure(5),18,18,18)
set(gca,'visible','off')

   subplot(1,2,2)
xlim([-50 50]);ylim([-60 60]);
axis square
Set_fig_YS(figure(5),18,18,18)
set(gca,'visible','off')

