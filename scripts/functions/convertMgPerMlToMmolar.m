function res = convertMgPerMlToMmolar(gpl,mw)
%UNTITLED3 Summary of this function goes here
%       gpl: milligrams/milliliter or grams/liter
%       mw: g/mole, or mg/millimole
%       ml: milliliters
%       res = millimoles/liter

milligrams = gpl;
molar = milligrams/mw;
res = molar * 1000;
end

