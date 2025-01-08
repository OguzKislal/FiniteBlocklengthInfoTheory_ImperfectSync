function [X,e] = gen_X_e_data(d,ts,tp,xm_p,dmax,SNR_lin)
% For a given delay, find X and e matrices from notes

% d is the delay, ts is the sampling period, tp is the pulse length 
% for the pilot sequence, xm_p is the pilot symbol sequence
% dmax is maximum delay

eps = rem(d,ts);
q = floor((d+1e-10)/ts); %1e-10 is added for numerical stabiltity.
N = tp/ts;
np = length(xm_p);

% total length of received sequence
M = dmax./ts +(np*tp)./ts +1;

xpN = ones(N,np);

for i = 1:np
    xpN(:,i) = xpN(:,i)*xm_p(i);
end

xpN = xpN(:);

%% Create upsampled data symbols for before and after pilot symbols

% length of padding at beginning for X1 and X2
q1 = q;
q2 = q + 1;

% length of padding at end for X1 and X2
r1 = M-q1-length(xpN);
r2 = M-q2-length(xpN);

num_symbols_start = ceil(q2/N);
num_symbols_end = ceil(r1/N);

start_symbols = zeros(num_symbols_start,1) ; 
% start_symbols = sqrt(SNR_lin).*( 2*(rand(num_symbols_start,1)>.5) - 1);
end_symbols = sqrt(SNR_lin).*(2*(rand(num_symbols_end,1)>.5) - 1);

% pad1_start = 2*(rand(q1,1)>.5) - 1;
% pad2_start = [pad1_start; 2*(rand(1,1)>.5) - 1];
% 
% pad1_end = 2*(rand(M-q1-length(xpN),1)>.5) - 1;
% pad2_end = pad1_end(1:end-1);

start_up = upsample(start_symbols,N);
end_up = upsample(end_symbols,N);

pad1_start = start_up(end-(q1-1):end);
pad2_start = start_up(end-(q2-1):end);

pad1_end = end_up(1:r1);
pad2_end = end_up(1:r2);

X1 = [pad1_start;xpN;pad1_end];
X2 = [pad2_start;xpN;pad2_end];

% X1 = [zeros(q1,1);xpN;zeros(M-q1-length(xpN),1)];
% X2 = [zeros(q2,1);xpN;zeros(M-q2-length(xpN),1)];

X = [X1 X2];
e = [1 - eps./ts ; eps./ts];
end
