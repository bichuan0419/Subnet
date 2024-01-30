function lr = likelihood_tuning(W_temp, CID_temp, r)
% this calculate the likelihood function of parameter tuning
% Chuan
%%
non_sing_loc = find(CID_temp~=1);
CID_len = length(non_sing_loc);
lr_vec = zeros(CID_len,1);
Clist_temp = 1:length(W_temp);
W_vec = squareform(W_temp);
edgeLen = length(W_vec);
global_dens = sum(squareform(W_temp)>r)/edgeLen;
for i = 1:CID_len
    sublist = getSubvector(Clist_temp, CID_temp,non_sing_loc(i));
    vec_idx = node2vec(sublist, sublist, length(Clist_temp));
    rest_idx = setdiff(1:edgeLen, vec_idx);
    cov_related_edges = W_vec(vec_idx);
    non_cov_related_edges = W_vec(rest_idx);
    pi_1 = sum(cov_related_edges > r)/(nchoosek(CID_temp(non_sing_loc(i)),2));
    pi_0 = sum(non_cov_related_edges>r)/(edgeLen - nchoosek(CID_temp(non_sing_loc(1)),2));
    e_beta1 = cov_related_edges > r;
    e_beta2 = non_cov_related_edges>r;
    term1 = sum(e_beta1*log(pi_1/global_dens+eps) +...
                (ones(size(e_beta1))- e_beta1)*log((1-pi_1)/(1-global_dens)+eps));
    term2 = sum(e_beta2(1,:)*log(pi_0/global_dens+eps) +...
                (ones(size(e_beta2))- e_beta2)*log((1-pi_0)/(1-global_dens)+eps));
    lr_vec(i) = term1 + term2;
end
lr = sum(lr_vec);
end