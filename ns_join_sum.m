function cellsum = ns_join_sum(celllist)

cellsum = cell2mat(celllist(1));
for i=2:length(celllist)
  cellsum = cellsum + cell2mat(celllist(i));
end

