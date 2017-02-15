function [avg_alpha, avg_theta, avg_phi] = lda(docs,K)






max_iters = 100;

collect_burn_in = 80;

collect_interval = 10;

hyper_burnin_in = 10;

alpha_sampling_interval = 2;

is_sym_alpha = false;

is_display = true;



N = max(docs(:,1));


V = max(docs(:,2));



n_doc_topic = zeros(N,K);

n_topic_word = zeros(K,V);

n_doc_dot = zeros(1,N);


n_topic_dot = zeros(1,K);

ii = docs(:,1);
vv = docs(:,2);
mm = docs(:,3);

num_words = sum(mm);

L = length(ii);

z = zeros(L, max(mm));

for l = 1:L
    i = ii(l);
    
    v = vv(l);
    
    m = mm(l);
    
    k = randi(K, 1, m);
    
    
    
    
    for ww = 1:m
        
        n_doc_topic(i,k(ww)) = n_doc_topic(i,k(ww)) + 1;
        n_doc_dot(i) = n_doc_dot(i) + 1;
        n_topic_word(k(ww),v) = n_topic_word(k(ww),v) + 1;
        n_topic_dot(k(ww)) =  n_topic_dot(k(ww)) + 1;
        
        z(l, ww) = k(ww);
        
    end
    
    
end



alpha = 0.1 * ones(1, K);

gamma = 0.01 * ones(1, V);

sum_alpha = sum(alpha);

sum_gamma = sum(gamma);

avg_alpha = 0;


avg_phi = 0;

avg_theta = 0;


avg_count = 0;

for r = 1:max_iters
    for l = 1:L
        i = ii(l);
        
        v = vv(l);
        
        m = mm(l);
        
        k = z(l,:);
        
        
        for ww = 1:m
            
            k_ww = k(ww);
            
            n_doc_topic(i,k_ww) = n_doc_topic(i,k_ww) - 1;
            n_doc_dot(i) = n_doc_dot(i) - 1;
            n_topic_word(k_ww,v) = n_topic_word(k_ww,v) - 1;
            n_topic_dot(k_ww) =  n_topic_dot(k_ww) - 1;
            
            
            
            p_left = (alpha + n_doc_topic(i,:));
            p_right = (gamma(v) + n_topic_word(:,v)) ./ (sum_gamma + n_topic_dot');
            p = (p_left .* p_right');
            
            sum_cum = cumsum(p(:));
            
            new_k_m = find(sum_cum > rand() * sum_cum(end),1);
            
            
            
            
            n_doc_topic(i,new_k_m) = n_doc_topic(i,new_k_m) + 1;
            n_doc_dot(i) = n_doc_dot(i) + 1;
            n_topic_word(new_k_m,v) = n_topic_word(new_k_m,v) + 1;
            n_topic_dot(new_k_m) =  n_topic_dot(new_k_m) + 1;
            
            z(l,ww) = new_k_m;
            
            
        end
        
        
    end
    
    
    
    
    if mod(r,10) == 0 && is_display
        temp_theta = (alpha + n_doc_topic) ./ (sum_alpha + n_doc_dot');
        temp_phi = (gamma + n_topic_word) ./ (sum_gamma + n_topic_dot');
        prob = sum(temp_theta(ii,:) .* temp_phi(:,vv)',2);
        pp = log(prob) .* mm;
        
        pp = exp(-sum(pp) ./ num_words);
        
        
        disp(sprintf('%d:%f\n',r,pp));
        
    end
    
    
    if alpha_sampling_interval > 0 && r >= hyper_burnin_in && mod(r,alpha_sampling_interval) == 0
        
        [alpha,sum_alpha] = sample_dir_concentration(n_doc_topic, alpha, is_sym_alpha);
        
    end
    
    
    
    if mod(r,collect_interval) == 0 && r >= collect_burn_in
        
        avg_alpha = avg_alpha + alpha;
        
        
        avg_theta = avg_theta + (alpha + n_doc_topic) ./ (sum_alpha + n_doc_dot');
        
        avg_phi = avg_phi + (gamma + n_topic_word) ./ (sum_gamma + n_topic_dot');
        
        avg_count = avg_count + 1;
        
    end
    
end

if avg_count > 0
    avg_alpha = avg_alpha ./ avg_count;
    avg_theta = avg_theta ./ avg_count;
    avg_phi = avg_phi ./ avg_count;
else
    avg_alpha = alpha;
    
    avg_theta = (alpha + n_doc_topic) ./ (sum_alpha + n_doc_dot');
    
    avg_phi = (gamma + n_topic_word) ./ (sum_gamma + n_topic_dot');
    
end

end






