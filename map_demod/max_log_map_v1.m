function Le = max_log_map_v1(y, var, La, symbols, labels)

    sym_len = length(y);
    y = reshape(y, 1, 1, sym_len);
    M = length(symbols);
    K = log2(M);
    La = reshape(La, 1, K, sym_len);
    L = zeros(K, sym_len);
    bin_labels = dec2bin(labels, K) - '0';

    for i = 1:K
        index_0 = bin_labels(:, i)==0;
        index_1 = bin_labels(:, i)==1;
        tmp_0 = max(-abs(y-symbols(index_0)).^2./var + ...
                    sum((1-bin_labels(index_0, :)).*La, 2), [], 1);
        tmp_1 = max(-abs(y-symbols(index_1)).^2./var + ...
                    sum((1-bin_labels(index_1, :)).*La, 2), [], 1);
        L(i, :) = tmp_0(:)-tmp_1(:);
    end
    Le = L(:)-La(:);
end
