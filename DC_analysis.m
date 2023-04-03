DC_1 = intersect(find(all_cluster==9),find(age_num==1));
DC_2 = intersect(find(all_cluster==11),find(age_num==1));
DC_all = [DC_1; DC_2];
all_cluster_temp = all_cluster;
all_cluster_temp(DC_all) = 0;
DC_not =intersect(find(all_cluster_temp),find(age_num==1));

Xcr1 = find(strcmp('Xcr1',name_joint)); % cDC1
Cd209a = find(strcmp('Cd209a',name_joint)); % cDC2
Sirpa = find(strcmp('Sirpa',name_joint)); % cDC2
Itgax = find(strcmp('Itgax',name_joint)); % cDC2
Itgam = find(strcmp('Itgam',name_joint)); % cDC2
Clec9a = find(strcmp('Clec9a',name_joint)); % cDC2
Cd24a = find(strcmp('Cd24a',name_joint)); % cDC2
Cd4 = find(strcmp('Cd4',name_joint)); % cDC2
Cd8a = find(strcmp('Cd8a',name_joint)); % cDC2
Cd14 = find(strcmp('Cd14',name_joint)); % cDC2
Cd163 = find(strcmp('Cd163',name_joint)); % cDC2
Clec10a = find(strcmp('Clec10a',name_joint)); % cDC2
Cx3cr1 = find(strcmp('Cx3cr1',name_joint)); % cDC2
Btla = find(strcmp('Btla',name_joint)); % cDC2
Cadm1 = find(strcmp('Cadm1',name_joint)); % cDC2


Gzmk = find(strcmp('Gzmk',name_joint))



DC_genes = [Itgax Xcr1 Clec9a Cd24a Cd8a Btla Cadm1 Cd209a Sirpa  Itgam   Cd4  Cd14 Cd163 Clec10a Cx3cr1  ];
DC_genes_red = [Itgax Xcr1 Clec9a Cd24a Btla Cd209a Sirpa  Itgam Clec10a ];

DC_1_mat = data_both_expressed{1}(DC_genes_red,DC_1);
DC_2_mat = data_both_expressed{1}(DC_genes_red,DC_2);
DC_not_mat = data_both_expressed{1}(DC_genes_red,DC_not);

DC_1_avg = mean(DC_1_mat,2);
DC_2_avg = mean(DC_2_mat,2);
DC_not_avg = mean(DC_not_mat,2);

%%

violin_markers = DC_genes_red;
marker_name = {'CD11c','Xcr1','Clec9a','CD24a','CD8a','Btla','CADM1','CD209','SIRPA','CD11b','CD4','CD14','CD163','Clec10a','CX3CR1'};
marker_name_red = {'CD11c','Xcr1','Clec9a','CD24a','Btla','CD209','SIRPA','CD11b','Clec10a'};

DC_mat_all = {DC_1_mat;DC_2_mat;DC_not_mat};

cluster_number = max(all_cluster);
num = 1:cluster_number;

for cnt_gene = 1:length(violin_markers)
    figure(1);
       subplot(9,1,cnt_gene); hold on
    
    gene_max = max(log(data_both_expressed{1}(violin_markers(cnt_gene),:)+1));
    gene_min = min(log(data_both_expressed{1}(violin_markers(cnt_gene),:)+1));
    k=1;
    for cnt_cluster = 1:3
        
        data = DC_mat_all{cnt_cluster}(cnt_gene,:);
        expressing = sum(data>0)/length(data);
        data_plot = ((log(data'+1)-gene_min)/(gene_max-gene_min))*expressing;
        violinplot_ys_lite(data_plot,num(k),0,'width',0.35,'ShowData',false,'ViolinColor',color(cnt_cluster,:),'MedianColor',color(cnt_cluster,:),'ShowMean',true);
        k=k+1;
        max_data_plot(cnt_gene,cnt_cluster) = max(data_plot);
    end
    
    top_value(cnt_gene) = round(max(max_data_plot(cnt_gene,:)),2);
%     set(gca,'xtick',[1:3],'xticklabels',[{'DC a','DC b', 'not DC'}])
set(gca,'xtick',[1:3],'xticklabels',[])
    set(gca,'ytick',[0 top_value(cnt_gene)],'yticklabels',[]);
    xlim([0 3.5])
%     ylabel(marker_name{cnt_gene})
end

%%

fraction_DC1_y = 100*(length(DC_1)/size(data_both_expressed{1},2))
fraction_DC2_y = 100*(length(DC_2)/size(data_both_expressed{1},2))

DC_1_o = intersect(find(all_cluster==9),find(age_num==2));
DC_2_o = intersect(find(all_cluster==11),find(age_num==2));

fraction_DC1_o = 100*(length(DC_1_o)/size(data_both_expressed{2},2))
fraction_DC2_o = 100*(length(DC_2_o)/size(data_both_expressed{2},2))

%% Spleen

spleen_CD3 = [32.12 35.64 34.52; 24.19 26.09 26.47];
cats = ({'Young','Old'});

figure(2);
violins_sp = violinplot(spleen_CD3', cats,'ViolinColor',[1 .7 .4],'Width',0.1,'ShowMean',true);
set(gca,'xtick',[1:2],'xticklabels',[{'Young','Old'}])
