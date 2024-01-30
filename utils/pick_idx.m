function [Z,Z_cell] = pick_idx(cluster_size, N)
% given a vector of cluster sizes, return the corresponding indices the
% vectorized edge list.


Z_cell = cell(length(cluster_size),1);

for i =1:length(cluster_size)
    temp = zeros(N);
    if i==1
        temp(1:cluster_size(i),1:cluster_size(i)) = 1;
    else
        temp(1+sum(cluster_size(1:i-1)):sum(cluster_size(1:i)), ...
            1+sum(cluster_size(1:i-1)):sum(cluster_size(1:i))) = 1;
    end
    temp = temp - diag(diag(temp));
    Z_cell{i} = find(squareform(temp));
end
Z = [];
for i = 1:length(cluster_size)
    Z = [Z, Z_cell{i}];
end

end