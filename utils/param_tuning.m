function [lambda_out, cut_out] = param_tuning(Wp, method, prctile_vec)
    % parameter tuning function
    % Input: 
    % Wp: Original shuffled (observed) adjacency matrix (-log(p-values))
    % lb: lower bound for lambda_0
    % ub: upper bound for lambda_0
    % By default, lb and ub are set to be lb = 0.5, ub = 0.7
    % and by default, the prior g(r) is assumed to be 5 distinct points
    
    %% User-defined parameters for parameter tuning
    num_in = nargin;
    if num_in == 1
        method = 'spectral';
        prctile_vec = [95,96,97,98,99];
    elseif num_in == 2
        prctile_vec = [95,96,97,98,99];
    end
    
    lam_vec = linspace(0.5,0.7,10);
    % For convenience, we 5 distinct values to approximate a uniformly distributed
    % prior distirbution
    % locate 90, 92 94 96 98th percentile for cut-off
    nlogp = squareform(Wp);
    gr = prctile(nlogp, prctile_vec);
    nodeLen = length(Wp);
    edgeLen =nodeLen*(nodeLen-1)/2;
    term_1_vec = zeros(length(lam_vec),length(gr));
    term_2_vec = zeros(length(lam_vec),length(gr));
    %% selection of lambda0 begins here
    for i = 1:length(lam_vec)
        
        lambda = lam_vec(i);
        % to speed up, we use k_means_iter = 1
        kmeans_iter = 1;

        for j = 1:length(gr)
            r = gr(j);
            
            % first determine the global edge density
            global_dens = sum(squareform(Wp)>r)/edgeLen;
            if strcmp(method, 'spectral')
                [CID_temp,W_temp, ~]=SICERS_final(Wp,r,lambda, kmeans_iter);
            elseif strcmp(method,'greedy')
                [W_temp, ~, CID_temp] = greedy_final(Wp, r, lambda);
            end
            term_sum(i,j) = likelihood_tuning(W_temp, CID_temp, r);
        end
        

    end
    %% choose lambda
%     term_sum = term_1_vec + term_2_vec;
    lr_lambda = gr*term_sum';
    [~,lambda_idx] = max(lr_lambda);
    lambda_idx = lambda_idx(1);
    lambda_out = lam_vec(lambda_idx);
    %% choose cut
    lr_cut = term_sum(lambda_idx,:);
    [~,cut_idx] = max(lr_cut);
    cut_idx = cut_idx(1);
    cut_out = gr(cut_idx);
end