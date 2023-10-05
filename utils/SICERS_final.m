function [CIDA,W_SICERS, Clist]=SICERS_final(W,r,lambda, kmeans_iter)
%%%% This function is for parsimonious detector


%%% Inputs:
%%%%%  W:   a n by n matrix of t test -log(p)-values from the raw data,
%%%%%  where n denotes the number of nodes.

%%%%%  p0:      A threshold on the -log(p)-values. We only do the clustering on the
%%%%%           significant edges.

%%% Outputs:
%%%%%  CIDA:    the cluster sizes of every cluster in a power descending
%%%%%           order. i.e. CID(1) will be the cluster size of the most
%%%%%           concentrated cluster
%%%%%  W_SICERS: Output reordered matrix
%%%%%  Clist:   the reordered node index


%%% Toy example:
%pattern one: one 20x20 corner
% % % % % Z90_2=eye(90);
% % % % % Z90_1=Z90_2;
% % % % % cluster_size = 20;
% % % % % Z90_1(1:cluster_size,1:cluster_size)=ones(cluster_size);
% % % % % figure;imagesc(Z90_1-Z90_2)
% % % % % a1=find(squareform(Z90_1-Z90_2)==1);
% % % % % a0=find(squareform(Z90_1-Z90_2)==0);
% % % % % Wt=zeros(4005,30);
% % % % % Wt(a0,:)=sqrt(5).*randn(4005-(cluster_size^2-cluster_size)/2,30);
% % % % % Wt(a1,:)=sqrt(5).*randn((cluster_size^2-cluster_size)/2,30)+0.8;
% % % % % Wt2=randn(4005,30);
% % % % % [h1,p1]=ttest(Wt',Wt2');
% % % % % ar=randperm(90);
% % % % % W1=squareform(-log(p1));
% % % % % W1=W1(ar,ar);
% % % % % nlogp=squareform(W1);
% % % % % 
% % % % % figure;imagesc(W1);
% % % % % [lambda_out, cut_out] = param_tuning(W1);
% % % % % [CIDA,W_SICERS, Clist]=SICERS_final(W1,cut_out,lambda_out, 5);
% % % % % figure;imagesc(W_SICERS)

%% Preprocessing of data
W1=W;
W(W1<r)=0;%Threshold on the p-values
z1=find(sum(W)>0); %Exclude the isolated nodes
W=W(z1,z1);
degs = sum(W, 2);
D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
L=D-W;%Laplacian matrix

%% Determine the number of clusters K
[U,d] = eig(L);
[d,ind] = sort(diag(abs(d)));
U = U(:,ind);
Mk=[];
Qual=[];

for m=1:kmeans_iter
    Prp_net=zeros(size(L,1),1);
    for K=1:5:size(L,1)
        %    C=kmeans(U(:,1:K),K,'Options',options,'Replicates',kmeans_iter);
        C=kmeans(U(:,1:K),K,'Replicates',kmeans_iter);
        
        indx=[]; %indx
        A_net=[];% chi2 in the net
        net_V=[];% size of each cluster
        C_net=[];
        for k=1:K
            indx=[indx;find(C==k)];
            net_V(k)=length(find(C==k));
            WC=W(find(C==k), find(C==k));
            C_net(k) = sum(WC(find(WC>0)))/2;
            A_net(k)=(net_V(k)*(net_V(k)-1))/2;
        end
        % (3.13) in supp
        Prp_net(K)=(sum(C_net)^(1-lambda))*(sum(C_net)/sum(A_net))^lambda;
    end
    
    Kbest = find(Prp_net == max(Prp_net));
    Kbest = Kbest(1);%In case several k's give the same Prp_net value
    Mk(m,:)=[Kbest max(Prp_net)];
    
    Qual(:,m)=Prp_net;
end

Kfinal=Mk(find(Mk(:,2)==max(Mk(:,2))),1);
Kfinal=Kfinal(1);


%% Find the cluster ID for each of the nodes
% [U,~] = eig(L);
C=kmeans(U(:,1: Kfinal),Kfinal,'Replicates',kmeans_iter);
indx=[];
A_net=[];
net_V=[];
C_net=[];
for k=1:Kfinal
    indx=[indx;find(C==k)];
    net_V(k)=length(find(C==k));
    WC=W(find(C==k), find(C==k));
    C_net(k) = sum(WC(find(WC>0)))/2;
    A_net(k)=(net_V(k)*(net_V(k)-1))/2;
end


diagscore=(C_net).^(1-lambda).*(C_net./A_net).*lambda;
diagscore(isnan(diagscore))=0;
[diagscore_sort,diagscore_sortID]=sort(diagscore,'descend');

CIDA = [];
inx_imporance=[];
for i=1:Kfinal
    inx_imporance=[  inx_imporance; find(C==diagscore_sortID(i))];
    CIDA = [CIDA; length(find(C==diagscore_sortID(i)))];
end

Cindx = 1:size(W1,1);
Cindx(z1)=C;
Cindx(setdiff(1:size(W1,1),z1))=-1;
CID=diagscore_sortID;
Clist = z1(inx_imporance);
Clist = [Clist setdiff(1:size(W1,1),z1)];
CIDA = [CIDA;length(setdiff(1:size(W1,1),z1))];
W_SICERS = W1(Clist,Clist);
end
%