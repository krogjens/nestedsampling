function ns_print(results,models,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prints out a summary of 'results' to the text file ['misc.data_id',misc.nssummary]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = [misc.data_id,misc.nssummary];
fid = fopen(fname,'w');

% Find best model (choose one)
for i = 1:length(models)
    evi(i) = results(i).Z_norm;
end
[~,best] = max(evi);
fprintf(fid,'Highest evidence for Model %i (probability %.2f). \n',best,results(best).Z_norm);

for i=1:length(models)
    fprintf(fid,'\n');
    fprintf(fid,'Model %i (probability %.3f, %i variable parameter',i,evi(i),models(i).options.lengthu);
    if models(i).options.lengthu==1
      fprintf(fid,'):\n');
    else
      fprintf(fid,'s):\n');
    end
    fprintf(fid,'Estimated value for log10-evidence: %.3f +- %.3f.\n',results(i).logZ(1)/log(10),results(i).logZ_error/log(10));
    fprintf(fid,'P-values [ H(replicated obs) > H(obs)]: \n');
    for j = 1:length(results(i).prob)
      fprintf(fid,'  n = %i:  %.2f. \n',j,results(i).prob(j));
    end 
    fprintf(fid,misc.titles{1});
    for j=1:length(misc.percentiles_at)
      fprintf(fid,'% .2f	',misc.percentiles_at(j));
    end
    fprintf(fid,[misc.titles{2},'	',misc.titles{3},'\n']);
    for j=1:length(models(i).labels)
        fprintf(fid,misc.labels(models(i).labels(j),:));
        for k=1:length(misc.percentiles_at)
          fprintf(fid,'% .3f	',results(i).percentiles(j,k));
        end
        fprintf(fid,'% .3f	% .3f	+/- % .3f\n',results(i).maxLpar(j),results(i).param_mean(j),results(i).param_stddev(j));
    end
    fprintf(fid,'log10-Maximal likelihood: % .3f\n',results(i).samples(end).logl/log(10));
end
fclose(fid);

