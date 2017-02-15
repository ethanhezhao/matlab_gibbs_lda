function cur_sum_value = learn_dir_sym_concentration(count_mat, cur_sum_value)
    
    num_dim = size(count_mat,2);
    
    [c_hist, v_hist]=hist(count_mat(:),unique(count_mat(:)));
    
    
    sum_over_col = sum(count_mat,2);
    
    [c_len,v_len]=hist(sum_over_col(:),unique(sum_over_col(:)));
    
    
    for iter = 1:200
        cur_parameter = cur_sum_value / num_dim;
        cur_digamma = 0;
        numerator = 0;
        
        start_idx = (v_hist(1) == 0) + 1; 
        for idx = start_idx:v_hist(end)
            cur_digamma = cur_digamma + 1.0 / (cur_parameter + idx - start_idx);
            count = find(v_hist == idx);
            if ~isempty(count)
                numerator = numerator + c_hist(count) * cur_digamma;
            end
        end
        
        cur_digamma = 0;
        
        denominator = 0;
        
        pre_len = 0;
        
        cached_digamma = psi(cur_sum_value);
        
        for idx = 1:length(c_len)
            len = v_len(idx);
            if len - pre_len > 20
                cur_digamma = psi(cur_sum_value + len) - cached_digamma;
            else
                for idxx = pre_len:len - 1
                    cur_digamma = cur_digamma + 1.0/ (cur_sum_value + idxx);
                end
            end
            
            denominator = denominator + cur_digamma * c_len(idx);
        end
        
        cur_sum_value = cur_parameter * numerator / denominator;
        
        assert(~isnan(cur_sum_value));
        
    end
end