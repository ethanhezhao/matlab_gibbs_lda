function [concentration,sum_con] = sample_dir_concentration(n_doc_topic, concentration , is_sym)
    sum_con = sum(concentration);
    K = length(concentration);
    if is_sym
        sum_con = learn_dir_sym_concentration(n_doc_topic, sum_con);
        concentration = sum_con / K * ones(1,K);
    else
        concentration = learn_dir_asym_concentration(n_doc_topic, concentration);
        sum_con = sum(concentration);
    end

end