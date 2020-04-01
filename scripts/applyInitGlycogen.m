function layout = applyInitGlycogen(layout,v)
%ApplyInitGlycogen(CometsLayout,struct) Sets an amount of glycogen relative
%to the initial biomass of the YeastGEM model, as set by parameters in v
%v must include fields: gly_molarmass, initgly_pct 
%v may include fields: initialpop, diff_gly

%if the fields don't exist, exit
if ~isfield(v,'gly_molarmass') || ~isfield(v,'initgly_pct')
    warning('Tried to invoke applyInitGlycogen, but did not provide v.gly_molarmass and v.initgly_pct');
    return
end
%if there's no glycogen, exit
if v.initgly_pct == 0
    return
end

%if biomass is already set, look it up. Otherwise use the parameter from v
x = ceil(layout.xdim / 2);
y = ceil(layout.ydim / 2);
p = 1e-6;
if isfield(v,'initialpop')
    p = v.initialpop;
end
if any(layout.initial_pop)
    %scan the initial pop field to find all colonies
    x = [];
    y = [];
    p = [];
    [~,xmax,ymax] = size(layout.initial_pop);
    for i = 1:xmax
        for j = 1:ymax
            if layout.initial_pop(1,i,j) > 0
                x = [x i];
                y = [y j];
                p = [p layout.initial_pop(1,i,j)];
            end
        end
    end
end

mass = v.gly_molarmass;
pct = v.initgly_pct; 
%convert pct to fractional weight, not percent
if pct > 1
    pct = pct / 100;
end

%get the glycogen concentration
millimolesglycogen = p * pct * 1000 / mass;

for i = 1:length(x)
    layout = setInitialMediaInCell(layout,x(i),y(i),'glycogen[c]',millimolesglycogen(i));
end

diffrate = 0;
if isfield(v,'diff_gly')
    diffrate = v.diff_gly;
end
layout = setDiffusion(layout,'glycogen[c]',diffrate);

