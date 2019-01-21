
function copyto = copyfield(copyto, copyfrom)
% copy fields from struct copyfrom to struct copyto
% WTJ, 20181017
if length(copyfrom) < length(copyto)
    error('copyfrom shorter than copyto!');
end
fns = fieldnames(copyfrom);
for jj = 1:length(fns)
    flds = [copyfrom.(fns{jj})];
    flds = num2cell(flds(1:length(copyto)));
    [copyto.(fns{jj})] = flds{:};
end


end
