function [v] = CR3BP_Stability(phi)
% Stability indices
% This function doesn't work 100% of the time and shouldn't be trusted.

eVal = eig(phi);

eValcomp = eVal(imag(eVal)~=0);

if isempty(eValcomp)
    eValcomp(1) = min(eVal);
    eValcomp(2) = 1/min(eVal);
end
    
eValreci = max(abs(eVal));

v1 = abs(eValreci + 1/eValreci)/2;
v2 = abs(eValcomp(1) + eValcomp(2))/2;

v = [v1 v2];

end

