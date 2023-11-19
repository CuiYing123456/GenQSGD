function [ set_out ] = DataTrans( set_in, row_num, col_num )
% data transportation for display and reshape format
set_size = size(set_in, 2);

for i = 1 : set_size
    X_tmp = set_in(:, i);
    X_tmp = reshape(X_tmp, row_num, col_num);
    X_tmp = X_tmp';
    X_tmp = reshape(X_tmp, row_num * col_num, 1);
    set_out(:, i) = X_tmp;
end
end

