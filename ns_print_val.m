function txt = ns_print_val(val,len)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prints a value with a format depending on its magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if val==0
  txt=' 0';
elseif val==1
  txt=' 1';
elseif abs(val)>=0.01 & abs(val)<100
  txt=sprintf('% .3f',val);
else
  txt=sprintf('% .1e',val);
end
for j=1:(len-length(txt))
  txt=[txt ' '];
end

