function [T_vec, score_max] = test_statistics(Wp, CID, cut_off)

W_vec  = squareform(Wp);

% initialize
T_vec = zeros(length(CID),1);


% define denominators
edge_num = length(W_vec);
edge_num2 = zeros(length(CID),1);
qhat = zeros(length(CID),1);
phat = sum(W_vec > cut_off)/edge_num;

for j = 1:length(CID)
    if j == 1
        curr_subnet = Wp(1:CID(1), 1:CID(1));
    else
        curr_subnet = Wp(sum(CID(1:j-1))+1:sum(CID(1:j)), sum(CID(1:j-1))+1:sum(CID(1:j)));
    end
    cur_subnet_vec = squareform(curr_subnet);
    edge_num2(j) =  length(cur_subnet_vec);
end

% define numerators
    
for j = 1:length(CID)
    if j == 1
        curr_subnet = Wp(1:CID(1), 1:CID(1));
    else
        curr_subnet = Wp(sum(CID(1:j-1))+1:sum(CID(1:j)), sum(CID(1:j-1))+1:sum(CID(1:j)));
    end
    if curr_subnet == 0
        curr_subnet = 1e-16*[0 1; 1 0];
    end
    cur_subnet_vec = squareform(curr_subnet);
    qhat(j) = sum(cur_subnet_vec > cut_off)/edge_num2(j);
end

for j = 1:length(CID)
    q_temp = qhat(j);
    cur_CID = CID(j);
    term1 = -(4/(q_temp - phat)^2 + 4/3/(q_temp - phat))^(-1);
    T_vec(j) = exp(term1 * cur_CID^2)*2*length(Wp);

end
score_max = max(T_vec);


end