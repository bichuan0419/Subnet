function vec_idx = node2vec(ClistA, ClistB, nodeLen)
    temp_zero = zeros(nodeLen);
    temp_zero(ClistA,ClistB) = 1;  
    temp_zero = temp_zero + temp_zero';
    temp_zero = temp_zero - temp_zero.*eye(nodeLen);
    temp_vec = squareform(temp_zero);
    vec_idx = find(temp_vec);


end