function enz = fpu2mmol(fpu,kcat,km)
%FPU2MMOL Convert the given amount of Filter Paper Units to mmols of enzyme
%   This is not 100% accurate because I'm using a linear approximation of
%   the integral for the catalysis rate (S/KM+S for the range of S from
%   50-48)

%The FPU line passes through 2mg glucose x the dilution at 1 FPU
%FPU = .37 / dilution
d = .37 / fpu;

%2 = enz * kcat * 3600 * 49/49+KM
enz_standard = 2 / (kcat * 3600 * 49 / (49 + km));

enz = enz_standard / d;
end

