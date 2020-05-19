function res = convertMgPerMlToMmoles(gpl,mw,ml)
%UNTITLED3 Summary of this function goes here
%       gpl: milligrams/milliliter or grams/liter
%       mw: g/mole, or mg/millimole
%       ml: milliliters
%       res = millimoles

milligrams = gpl * ml;
res = milligrams/mw;

end

