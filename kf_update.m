function [Xh, P] = kf_update(Xh_, P_, Z, Klmn)
% function [Xh P] = kalmanUpdate(Xh_, P_, Z, Klmn)
%
% H is Klmn.H, linear measurement model
% R is Klmn.R, noise covariance matrix

% get state vector
nx = size(Xh_,1);

% measurement covariance
S = Klmn.H*P_*Klmn.H' + Klmn.R;

% projected estimate
Zh = Klmn.H*Xh_;

% innovation
inov = Z-Zh;

% Klmn gain
W = P_*Klmn.H'/S;

[~, msgid]=lastwarn;
if strcmp(msgid, 'MATLAB:illConditionedMatrix') || ...
        strcmp(msgid, 'MATLAB:singularMatrix')
%     keyboard;
end 

% state update
Xh = Xh_ + W*inov;

% covariance update
P = (eye(nx)-W*Klmn.H)*P_;