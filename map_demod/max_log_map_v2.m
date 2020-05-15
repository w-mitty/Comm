function Le = max_log_map_v2(y, var, La, symbols, labels)

    sym_len = length(y);
    y = reshape(y, 1, 1, sym_len);
    M = length(symbols);
    K = log2(M);
    La = reshape(La, 1, K, sym_len);
    L = zeros(K, sym_len);
    bin_labels = dec2bin(labels, K) - '0';
    
    index_0 = zeros(M, K);
    index_1 = zeros(M, K);
    for i = 1:K
        index_0(:,i) = bin_labels(:, i)==0;
        index_1(:,i) = bin_labels(:, i)==1;
    end
    index_0 = logical(index_0);
    index_1 = logical(index_1);
    
    symbols = repmat(symbols,K,1);
    bin_labels = repmat(bin_labels,K,1);
    
    tmp_0 = max(-abs(y-reshape(symbols(index_0),M/2,1,1,K)).^2./var + ...
            sum((1-permute(reshape(bin_labels(index_0, :)',K,M/2,1,K),[2,1,3,4])).*La, 2), [], 1);
    tmp_1 = max(-abs(y-reshape(symbols(index_1),M/2,1,1,K)).^2./var + ...
                sum((1-permute(reshape(bin_labels(index_1, :)',K,M/2,1,K),[2,1,3,4])).*La, 2), [], 1);

    L = squeeze(tmp_0-tmp_1)';

    Le = L(:)-La(:);
    
end





