function Le = MAP_DEM(y,sigma,consteWithLabel,La)

    %notice: L = p0/p1
    
    L = length(y);
    M = size(consteWithLabel,1);
    K = log2(M);
    Le = zeros(size(La));
    for l = 1:L
        for k = 1:K

            p1 = 0;
            p0 = 0;
            index_1=consteWithLabel(:,k+1)==1;
            index_0=consteWithLabel(:,k+1)==0;
            subset_1=consteWithLabel(index_1,:);
            subset_0=consteWithLabel(index_0,:);

            for i = 1:M/2 
                
                p1 = p1+exp(-abs(subset_1(i,1)-y(l))^2/(2*sigma^2)+dot(La((l-1)*K+1:l*K),1-subset_1(i,2:end)));
                p0 = p0+exp(-abs(subset_0(i,1)-y(l))^2/(2*sigma^2)+dot(La((l-1)*K+1:l*K),1-subset_0(i,2:end)));

            end
            
            Le((l-1)*K+k) = log(p0/p1)-La((l-1)*K+k);

        end
    end
end
