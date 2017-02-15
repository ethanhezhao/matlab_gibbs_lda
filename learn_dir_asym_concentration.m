function cur_params = learn_dir_asym_concentration(count_mat, cur_params)

num_iter = 200;

shape = 1.00001;

scale = 1.0;

K = size(count_mat,2);

[c_hist, v_hist]=hist(count_mat,unique(count_mat));

sum_over_col = sum(count_mat,2);

[c_len,v_len]=hist(sum_over_col(:),unique(sum_over_col(:)));

for iter = 1:num_iter

    denominator = 0;
    
    cur_digamma = 0;
    
    sum_param = sum(cur_params);
    
    start_idx = (v_len(1) == 0) + 1; 

    for idx = start_idx:v_len(end)
        cur_digamma = cur_digamma + 1 ./ (sum_param + idx - start_idx);
        len = find(v_len == idx);
        if ~isempty(len)
            denominator = denominator + c_len(len) * cur_digamma;
        end
    end
    
    denominator = denominator - 1/scale;
    
    for k = 1:K
        
        nnz = find(c_hist(:,k));
        non_zero_limit = v_hist(nnz(end));
        
        cur_digamma = 0;
        
        old_params_k = cur_params(k);
        
        cur_params(k) = 0;
        
        start_idx = (v_hist(nnz(1)) == 0) + 1; 
        
        for idx = start_idx:non_zero_limit
            cur_digamma = cur_digamma + 1./(old_params_k + idx - start_idx);
            count = find(v_hist == idx);
            if ~isempty(count)
                cur_params(k) = cur_params(k) + c_hist(count,k) * cur_digamma;
            end
        end
        
        cur_params(k) = (old_params_k * cur_params(k) + shape) /denominator;
        
        assert(cur_params(k) >0);
    end

end