function [NM, AM, pdf, x] = PEDiri(a,Order,N,Nx)
%PEDIRI is designed to test analytic expressions of PE moments with
%Dirichlet priors by comparing their moments (up to the 4th) to that of 
%numerically-generated distributions. The output (Bayesian) probability
%density function is also computed. For prior distributions, set 'a' to be
%commensurate with prior knowledge (a = ones(1,D) for a uniform prior). 
%For posterior distributions add ordinal pattern counts from ORDENCODE to
%the prior 'a'.
%
%INPUTS
%
%a = Hyperparameters of the Dirichlet distribution.
%Order = Number of moments to compute (Max = 4).
%N = Number of samples to draw from the Dirichlet distribution.
%Nx = Number of elements in the domain of theoretically-calculated 
%distributions. 
%
%OUTPUTS
%
%NM = Numerically-computed moments of H.
%AM = Analytically-calculated moments of H.
%pdf = Theoretically calculated standard probability distributions using
%the computed moments. See below for a list of these distributions.
%x = The domain over which pdf is calculated.

%% Generate numerical samples

A = -1/log(length(a));
P = gamrnd(repmat(a,N,1),1,N,length(a));
P = P ./ repmat(sum(P,2),1,length(a));        %Each row of P is a sample
                                              %from the Dirichlet
                                              %distribution.
H = A*sum(P.*log(P),2);                       
LA = length(a);

%% Compute moments of H

NM(1) = mean(H);
NM(2) = var(H);
NM(3) = skewness(H,0);
NM(4) = kurtosis(H,0);

%% Calculate the mean of H

a0 = sum(a);
AM(1) = A*sum(a/a0.*dpsi(0,a+1,a0+1));

%% Calculate the variance of H

if Order > 1
    t2_1 = sum(a.*(a+1)/a0/(a0+1).*(dpsi(0,a+2,a0+2).^2 + dpsi(1,a+2,a0+2)));
    T2 = (a'*a)/a0/(a0+1).*(dpsi(0,a'+1,a0+2)*dpsi(0,a+1,a0+2) - psi(1,a0+2));
    t2_2 = sum(T2(:))-trace(T2);
    AM(2) = A^2*(t2_1 + t2_2) - AM(1)^2;
end

%% Calculate skewness of H

if Order > 2
    t3_1 = sum(a.*(a+1).*(a+2)/a0/(a0+1)/(a0+2).*(dpsi(0,a+3,a0+3).^3 +...
        3*dpsi(1,a+3,a0+3).*dpsi(0,a+3,a0+3) + dpsi(2,a+3,a0+3)));
    T3_1 = 3*(a'.*(a'+1))*a/a0/(a0+1)/(a0+2).*((dpsi(0,a'+2,a0+3).^2 +...
        dpsi(1,a'+2,a0+3))*dpsi(0,a+1,a0+3) -...
        repmat(2*dpsi(0,a'+2,a0+3)*psi(1,a0+3),1,LA) - ...
        psi(2,a0+3));
    t3_2 = sum(T3_1(:))-trace(T3_1);
    T3_2 = Outerprod3(a,a,a)/a0/(a0+1)/(a0+2).*(Outerprod3(dpsi(0,a+1,a0+3),dpsi(0,a+1,a0+3),dpsi(0,a+1,a0+3)) -...
        (repmat(dpsi(0,a'+1,a0+3),[1, LA, LA]) +...
        repmat(dpsi(0,a+1,a0+3),[LA, 1, LA]) +...
        repmat(permute(dpsi(0,a+1,a0+3),[1 3 2]),[LA, LA, 1]))*psi(1,a0+3) -...
        psi(2,a0+3));
    t3_3 = SumUniqueInd3(T3_2);
    AM(3) = (A^3*(t3_1 + t3_2 + t3_3) - 3*AM(1)*AM(2) - AM(1)^3)*AM(2)^(-3/2);
end

%% Calculate the kurtosis of H

if Order > 3
    A4_1 = a.*(a+1).*(a+2).*(a+3)/a0/(a0+1)/(a0+2)/(a0+3);
    T4_1 = dpsi(0,a+4,a0+4).^4 + 6*dpsi(0,a+4,a0+4).^2.*dpsi(1,a+4,a0+4) +...
            dpsi(0,a+4,a0+4).*dpsi(2,a+4,a0+4) +...
            3*dpsi(0,a+4,a0+4).*dpsi(2,a+4,a0+4) +...
            3*dpsi(1,a+4,a0+4).^2 +...
            dpsi(3,a+4,a0+4);
    t4_1 = sum(A4_1.*T4_1);

    A4_2 = (a'.*(a'+1))*(a.*(a+1))/a0/(a0+1)/(a0+2)/(a0+3);
    T4_2 = (dpsi(1,a'+2,a0+4)+dpsi(0,a'+2,a0+4).^2)*(dpsi(1,a+2,a0+4)+dpsi(0,a+2,a0+4).^2) -...
            4*(dpsi(0,a'+2,a0+4)*dpsi(0,a+2,a0+4))*psi(1,a0+4) + 2*psi(1,a0+4)^2 -...
            2*(repmat(dpsi(0,a'+2,a0+4),1,LA) + repmat(dpsi(0,a+2,a0+4),LA,1))*psi(2,a0+4) -...
            psi(3,a0+4);
    t4_2 = sum(A4_2(:).*T4_2(:)) - trace(A4_2.*T4_2);

    A4_3 = (a'.*(a'+1).*(a'+2))*a/a0/(a0+1)/(a0+2)/(a0+3);
    T4_3 = dpsi(0,a'+3,a0+4).^3*dpsi(0,a+1,a0+4) +... 
            3*(dpsi(0,a'+3,a0+4).*dpsi(1,a'+3,a0+4)*dpsi(0,a+1,a0+4) - repmat(dpsi(1,a'+3,a0+4),1,length(a))*psi(1,a0+4)) -...
            3*repmat(dpsi(0,a'+3,a0+4).^2,1,LA)*psi(1,a0+4) -...
            3*repmat(dpsi(0,a'+3,a0+4)*psi(2,a0+4),1,LA) +...
            dpsi(2,a'+3,a0+4)*dpsi(0,a+1,a0+4) - psi(3,a0+4);
    t4_3 = sum(A4_3(:).*T4_3(:)) - trace(A4_3.*T4_3);

    A4_4 = Outerprod3(a.*(a+1),a,a)/a0/(a0+1)/(a0+2)/(a0+3);
    T4_4 = Outerprod3(dpsi(1,a+2,a0+4)+dpsi(0,a+2,a0+4).^2,dpsi(0,a+1,a0+4),dpsi(0,a+1,a0+4)) -...
            psi(1,a0+4)*(repmat(dpsi(1,a'+2,a0+4)+dpsi(0,a'+2,a0+4).^2,[1, LA, LA]) +...
                    2*repmat(dpsi(0,a'+2,a0+4)*dpsi(0,a+1,a0+4),[1, 1, LA]) +...
                    2*repmat(permute(dpsi(0,a'+2,a0+4)*dpsi(0,a+1,a0+4),[1 3 2]),[1, LA, 1])) -...
            psi(2,a0+4)*(2*repmat(dpsi(0,a'+2,a0+4),[1, LA, LA]) +...
                    repmat(dpsi(0,a+1,a0+4),[LA, 1, LA]) +...
                    repmat(permute(dpsi(0,a'+1,a0+4),[3 2 1]),[LA, LA, 1])) +...
            2*psi(1,a0+4)^2 - psi(3,a0+4);
    t4_4 = SumUniqueInd3(A4_4.*T4_4);

    A4_5 = Outerprod4(a,a,a,a)/a0/(a0+1)/(a0+2)/(a0+3);
    T4_5 = Outerprod4(dpsi(0,a+1,a0+4),dpsi(0,a+1,a0+4),dpsi(0,a+1,a0+4),dpsi(0,a+1,a0+4)) -...
            psi(1,a0+4)*(repmat(dpsi(0,a'+1,a0+4)*dpsi(0,a+1,a0+4),[1, 1, LA, LA]) +...
                    repmat(permute(dpsi(0,a'+1,a0+4)*dpsi(0,a+1,a0+4),[1 3 2 4]),[1, LA, 1, LA]) +...
                    repmat(permute(dpsi(0,a'+1,a0+4)*dpsi(0,a+1,a0+4),[1 4 3 2]),[1, LA, LA, 1]) +...        
                    repmat(permute(dpsi(0,a'+1,a0+4)*dpsi(0,a+1,a0+4),[3 2 1 4]),[LA, 1, 1, LA]) +...
                    repmat(permute(dpsi(0,a'+1,a0+4)*dpsi(0,a+1,a0+4),[4 2 3 1]),[LA, 1, LA, 1]) +...
                    repmat(permute(dpsi(0,a'+1,a0+4)*dpsi(0,a+1,a0+4),[3 4 1 2]),[LA, LA, 1, 1])) -...
            psi(2,a0+4)*(repmat(dpsi(0,a'+1,a0+4),[1, LA, LA, LA]) +...
                    repmat(dpsi(0,a+1,a0+4),[LA, 1, LA, LA]) +...
                    repmat(permute(dpsi(0,a'+1,a0+4),[3 2 1 4]),[LA, LA, 1, LA]) +...
                    repmat(permute(dpsi(0,a'+1,a0+4),[4 2 3 1]),[LA, LA, LA, 1])) +...
            3*psi(1,a0+4)^2 - psi(3,a0+4);
    t4_5 = SumUniqueInd4(A4_5.*T4_5);

    AM(4) = (A.^4*(t4_1 + 3*t4_2 + 4*t4_3 + 6*t4_4 + t4_5) - 4*AM(1)*AM(3)*AM(2)^(3/2) -...
            6*AM(1)^2*AM(2) - AM(1)^4)/AM(2)^2;
end

%% Generate Probability Distribution

x = linspace(0,1,Nx);

a = AM(1)^2/AM(2)*(1-AM(1)) - AM(1);
b = a*(1/AM(1)-1);
pdf = x.^(a-1).*(1-x).^(b-1)*exp(gammaln(a+b)-gammaln(a)-gammaln(b));
pdf = pdf/sum(x(2)*pdf(not(isinf(pdf))));

end

function [x] = dpsi(n,a,b)
%DPSI takes the difference of two psi functions;

x = psi(n,a)-psi(n,b);

end

function [OP] = Outerprod3(a1, a2, a3)
%OUTERPROD3 computes the m-by-n-by-p outer product array of row vectors of
%length m, n and p. Each element of the array is a unique combination of 
%products taking one term from each vector. 

A = permute(a1'*a2,[3 2 1]);
OP = zeros(length(a3),length(a2),length(a1));

for i = 1:length(a1)
    OP(:,:,i) = a3'*A(:,:,i);
end

OP = permute(OP,[3 2 1]);

end

function [OP] = Outerprod4(a1, a2, a3, a4)
%OUTERPROD4 computes the m-by-n-by-p-by-q outer product array of row 
%vectors of length m, n, p and q. Each element of the array is a unique 
%combination of products taking one term from each vector. 

A = permute(Outerprod3(a1,a2,a3),[4 2 3 1]);
OP = zeros(length(a4),length(a2),length(a3),length(a1));

for i = 1:length(a3)
for j = 1:length(a1)
    OP(:,:,i,j) = a4'*A(1,:,i,j);
end
end

OP = permute(OP,[4 2 3 1]);

end

function [s] = SumUniqueInd3(A)
%SUMUNIQUEIND3 takes an m-by-m-by-m 3D array and adds all the elements with 
%indices that have no repeats.
%
%e.g. for a 3-by-3-by-3 matrix, the elements defined by the indices 123,
%213, 132, 231, 312, 321 would be included in the summation, but not 111, 
%112, 113 etc. since they contain repeating indices. 

S = zeros(1,size(A,1));

for i = 1:size(A,1)
    B = A(:,:,i);
    B(i,:) = [];
    B(:,i) = [];
    S(i) = sum(B(:))-trace(B);
end

s = sum(S);

end

function [s] = SumUniqueInd4(A)
%SUMUNIQUEIND4 takes an m-by-m-by-m-by-m 4D array and adds all the elements 
%that have no repeating indices.
%
%e.g. for a 3-by-3-by-3 matrix, the elements defined by the indices 123,
%213, 132, 231, 312, 321 would be included in the summation, but not 111, 
%112, 113 etc. since they contain repeating indices. 

S = zeros(size(A,1));

for i = 1:size(A,1)
for j = 1:size(A,1)
    B = A(:,:,i,j);
    B([i j],:) = [];
    B(:,[i j]) = [];
    S(i,j) = sum(B(:))-trace(B);
end
end

s = sum(S(:))-trace(S);

end