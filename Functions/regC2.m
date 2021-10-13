function  C = regC2(Csample)

I_N = size(Csample,1)/2;
V_N = I_N;

CsampleVV = Csample(1:V_N,1:V_N);
CsampleIV = Csample(V_N+1:V_N+I_N,1:V_N);
CsampleII = Csample(V_N+1:V_N+I_N, V_N+1:V_N+I_N);
C = zeros(I_N, V_N);
for i = 1:I_N    %the I value
    for j = 1:V_N   %the V value
        uIdx = setdiff(1:V_N, [i, j]);   %no V_i, V_j
        Sy1u = CsampleIV(i, uIdx);
        Sy2u = CsampleVV(j, uIdx);
        invUU = inv(CsampleVV(uIdx, uIdx));
        C(i,j) = CsampleIV(i,j) - Sy2u*invUU*Sy1u';
    end
end

