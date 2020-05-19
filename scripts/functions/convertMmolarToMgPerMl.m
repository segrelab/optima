function res = convertMmolarToMgPerMl(mmolar,mw)
%convertMmolarToMgPerMl converts the given number of millimoles of substance with molecular weight mw (g/mole) into milligrams per milliliter 
%   units:
%       mmolar: millimoles/milliliter
%       mw: g/mole, or mg/millimole
%       res: mg/mL

milligramsperliter = mw * (mmolar);
res = milligramsperliter / 1000;

% g/mole * millimoles / liters -> g/mole * moles/1000L -> 
%    g*mole / 1000moles*liters -> g/1000L = mg/L

end

