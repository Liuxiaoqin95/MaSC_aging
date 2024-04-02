%entropy analysis 

clc;clear;

stage=string(virgincount_metadata.stage);
clusterIfo.idxCluster={};
s1=find(stage=='s1');
s2=find(stage=='s2');
s3=find(stage=='s3');
s4=find(stage=='s4');
clusterIfo.idxCluster={s1 s2 s3 s4};

a=1:length(stage);
for k=1:length(stage)
    a(k)=str2num(strrep(stage(k),"s",""));
end

clusterIfo.identity=(a)';

%vp30count=vp30count(:,remaincell);
data=sparse(virgincount);
%genes
cells=(a)';

tmpIniData.data=data;
tmpIniData.genes=virgincount_gene;
tmpIniData.cells=cells;


minCells = 0; minGenes = 0; logNormalize = 1; filterRibo = 0; % default parameters for QC (see functions for details):
proData = preprocessing(tmpIniData,minCells, minGenes,logNormalize,filterRibo);
%   proData: a struct variable (store the data after QC) contains the same fields with iniData
%quick_construct = 0; tau = [];
%quick_construct = 1; tau = 0.4;

quick_construct = 1; tau =0.1;  %0.1..0.7

networkIfo = constructingNetwork(proData.data',quick_construct,tau); % a struct variable

[scE,scEcell] = estimatingscEnergy(proData.data,networkIfo);
% scE: energy matrix, the same dimension with data
% scEcell; scEnergy of each cell
% perform principal component analysis of the energy matrix scE
ydata = ECA(scE);% ydata : m x nPC, nPC-D coordindates from dimension reduction
%ydata = ydata(:,1:2); % In this demo, only the first two significant components are used. 
%ydata0=ydata;
%ydata = ydata0(:,1:10);

rootNode=1;
numCluster = length(unique(clusterIfo.identity));



alpha = 0.01; theta1 = 0.8; % default parameters (see functions for details):
lineageIfo = inferingLineage(scEcell,ydata,clusterIfo,rootNode,alpha,theta1);

Rscript = '"D:/software2/R-3.6.1/bin/Rscript"'; % for 64-bit windows
%Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
% User may also need to change the path of R script to execute inside inferingPseudotime.m if current directory is not ./scEpath-master/
pseudotimeIfo = inferingPseudotime(Rscript,ydata,lineageIfo,clusterIfo);


colorCell = distinguishable_colors(numCluster);% colors for each cluster
group = clusterIfo.identity; % m x 1 numerical vector, the cluster assignment for each cell
class_labels = strcat('C',cellstr(num2str([1:numCluster]'))); % text annotations of each cluster

if 1==2
% visualize cells on two-dimensional space
% true_labs = []; % it can be empty
true_labs = proData.cells;
true_labs = categorical(true_labs);
%true_labs = reordercats(true_labs,{'E14.5','E16.5','E18.5','Adult'});
marker_size = scEcell*50;  % the size of individual cell (dots), propotational to the scEnergy of each cell (by default)
fig_width = 600;fig_height = 250;
cluster_visualization(ydata, group,class_labels, true_labs, marker_size,colorCell,fig_width,fig_height)
end
node_size = grpstats(full(scEcell),group,'median')*30; % the size of tree node, propotational to the median scEnergy of each cluster
showLoops = 1; fig_width = 200;
lineage_visualization(lineageIfo,class_labels,node_size,colorCell,showLoops,fig_width,num2str(tau))

scEnergy_comparison_visualization(scEcell,clusterIfo,class_labels,colorCell,1,num2str(tau))

save("scEcell_1000gene_tau0.4.mat",'scEcell');
save("netid.mat","netid");
networkIfo_virgin=networkIfo;
save("networkIfo_virgin.mat","networkIfo_virgin");