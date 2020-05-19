function res = convertMmolesToMgPerMl(mmols,mw,ml)
%convertMmolesToMgPerMl converts the given number of millimoles of substance with molecular weight mw (g/mole) into milligrams per milliliter 
%   units:
%       mmols: millimoles
%       mw: g/mole, or mg/millimole
%       ml: milliliters
%       res: mg/mL

milligrams = mw * (mmols);
res = milligrams / ml;
end

