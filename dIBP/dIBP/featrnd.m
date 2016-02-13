function feature = featrnd(nobj, nfeat, p_thresh)
% Generate naive binary features without singular object
    feature = (rand(nobj, nfeat) > p_thresh);
    for i=1:nobj
        while (sum(feature(i, :)) == 0)
            feature(i, :) = (rand(1, nfeat) > p_thresh);
        end
    end
end