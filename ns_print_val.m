function txt = ns_print_val(val,len)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prints a value with a format depending on its magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(val,1)==0
  txt=sprintf('% -i',val);
elseif abs(val)>=0.01 & abs(val)<100
  if mod(val,0.01)==0
    txt=sprintf('% .2f ',val);
  else
    txt=sprintf('% .3f',val);
  end
else
  txt=sprintf('% .1e',val);
end
for j=1:(len-length(txt))
  txt=[txt ' '];
end

