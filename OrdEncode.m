function [symbols] = OrdEncode(y,D,tau,Overlap)
%==========================================================================
% DESCRIPTION
% Encodes a string of ordinal data into symbols [1-D!], where each symbol
% represents a distinct ordinal pattern.
%
% Version 1.0

% INPUTS
% y = Ordinal data. Must be a numeric data type.
% D = Ordinal pattern length (a.k.a. embedding dimension). This must be an
% integer in the range 2 <= D <= 8.
% tau = Partition number (a.k.a. embedding delay), encodes patterns based
% on every tau-th data point. Input ordinal data are separated into tau
% partitions, e.g. for tau = 2, odd-indexed and even-indexed data are
% placed in separate partitions. This must be a positive integer <
% floor(length(y)/D).
% Overlap = Boolean value denoting whether subsequent ordinal patterns are
% computed from overlapping data points. By default this is off, so that
% ordinal patterns are computed from e.g. indices 1-3, 4-6, 7-9 etc. With
% this option flagged, ordinal patterns are computed from indices 1-3, 2-4,
% 3-5, 4-6 etc. While overlapping subsets produce more symbols, they also
% introduce correlations into the symbol string.

% OUTPUTS
% symbols = Output symbol string (or cell array if the input tau was a 
% vector. Partitioned data are separated by a "0" (zero) symbol, e.g. for 
% tau = 2, the output will read cat(<encoded even indices> 0 <encoded odd 
% indices>).
%==========================================================================

%% Check inputs
assert(isnumeric(y), 'Ordinal data must be numeric')
assert(isvector(y), 'Ordinal data must be a vector')
assert(isnumeric(D), 'Dimension must be numeric')
assert(isscalar(D), 'Dimension must be a scalar')
assert(D == round(D) && D >= 2 && D <= 8, ...
    'Dimension must be an integer between 2 and 8 inclusive')
assert(isnumeric(tau), 'Embedding delay must be numeric')
assert(isvector(tau), 'Embedding delays must be contained in a vector')
assert(min(tau == round(tau)), 'Embedding delays must be integers')
assert(Overlap == 1 || Overlap == 0 || islogical(Overlap), ...
    'Overlap flag must be a logical value')

%% Begin Code
PERM = perms(1:D);
Dim = factorial(D);

%Preallocate matrices;
Ineq = zeros(Dim*(D-1),1);
A = zeros(1,sum(1:D-1));
B = zeros(1,sum(1:D-1));
symbols = cell(1,length(tau));

%% Construct/define the necessary inequality operations.
for i = 1:D-1                       
    Pairs((i-1)*Dim+1:i*Dim,:) = PERM(:,i:i+1);
end

if D/2 == ceil(D/2)                         %Check if D is even.
    for j = 1:(D/2-1)                       %For even D.
        A(1+D*(j-1):D*j) = 1:D;
        B(1+D*(j-1):D*j) = mod((1:D)+j-1,D)+1;
    end
    A(D*(D/2-1)+1:end) = 1:D/2;
    B(D*(D/2-1)+1:end) = (D/2+1):D;
else
    for j = 1:(D-1)/2                       %For odd D.
        A(1+D*(j-1):D*j) = mod(j*(1:D),D)+1;
        B(1+D*(j-1):D*j) = mod(j*((1:D)+1),D)+1;
    end
end

for k = 1:sum(1:D-1)
    Hit = (A(k) == Pairs(:,1)) & (B(k) == Pairs(:,2));
    Hit_Inv = (B(k) == Pairs(:,1)) & (A(k) == Pairs(:,2));
    Ineq(Hit) = k;              %Encode matching pairs as k.
    Ineq(Hit_Inv) = -1*k;       %Encode a match in the inverse order as -k.
end

RIneq = reshape(Ineq,Dim,[]);

%% Encode ordinal data into symbols.
for k = 1:length(tau)

    %Partition time series into a tau-rowed matrix;
    Y = zeros(D,length(y)+(1-D)*tau(k));           %Preallocate.
    
    for i = 1:D
        Y(i,:) = y((1+(i-1)*tau(k)):end-(D-i)*tau(k))';
    end
    if Overlap == 0
        Y(:,mod(0:size(Y,2)-1,D*tau(k)) > tau(k)-1) = [];
    end
    
    %Add random perturbations for the purpose of tie-breaking. 
    Y = Y + std(y)*1e-9*randn(size(Y));                 
    
    %Take all sum(1:D-1) logical operations;
    Logic_Y = Y(A,:) > Y(B,:);
    Logic_Ynot = not(Logic_Y);
    
    %Encode ordinal patterns based on the output of logical inequalities;
    OP = zeros(1,size(Y,2));
    for jj = 1:Dim
        Match = 1;
        for kk = 1:D-1
            if RIneq(jj,kk) > 0
                Match = and(Match,Logic_Y(RIneq(jj,kk),:));
            else
                Match = and(Match,Logic_Ynot(abs(RIneq(jj,kk)),:));
            end
        end
        OP(Match == 1) = jj;
    end
    
    %Organise symbols so that the symbol string from each partition is
    %separated by a "0" (zero).
    
    OPMat = [reshape([OP zeros(1,ceil(length(OP)/tau(k))*tau(k)-...
        length(OP))],tau(k),[])';zeros(1,tau(k))];
    OP = OPMat(:)';
    symbols{k} = OP;
    
end
