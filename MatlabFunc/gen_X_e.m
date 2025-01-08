function [X,e] = gen_X_e(d,ts,tp,xm_p,dmax)
% For a given delay, find X and e matrices from notes

% d is the delay, ts is the sampling period, tp is the pulse length 
% for the pilot sequence, xm_p is the pilot symbol sequence
% dmax is maximum delay

eps = rem(d,ts);
q = floor((d+1e-10)/ts); % 1e-10 is added for numerical stability.
N = tp/ts;
np = length(xm_p);

% total length of received sequence
M = dmax./ts +(np*tp)./ts +1;


xpN = ones(N,np);

for i = 1:np
    xpN(:,i) = xpN(:,i)*xm_p(i);
end

xpN = xpN(:);

q1 = q;
q2 = q + 1;


X1 = [zeros(q1,1);xpN;zeros(M-q1-length(xpN),1)];
X2 = [zeros(q2,1);xpN;zeros(M-q2-length(xpN),1)];

X = [X1 X2];
e = [1 - eps./ts ; eps./ts];
end