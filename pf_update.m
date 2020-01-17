function [p, wts] = pf_update(p, wts, Z, glfn)
% function [p wts] = pf_update(p, wts, Z, gh)
%
% p is d x N matrix with d dimensional state and N particles
% wts is 1 x N vector of weights for each particle
% Z is the measurement
% glfn is the measurement model

for jj=1:size(p,2)
      % ^^ each choice has different performance
%     wts(jj)=wts(jj).*glfn(Z, p(:,jj));  
    wts(jj)=glfn(Z, p(:,jj));
end

wts=wts/sum(wts);
neff=1/sum(wts.^2);
if (neff <=  size(p,2))
    p=p(:,resample(wts));
    wts=ones(1,numel(wts))/numel(wts);
end 