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
fprintf(fid,'Highest evidence for Model %i (probability %.3f). \n',best,results(best).Z_norm);

for i=1:length(models)
    fprintf(fid,'\n');
    fprintf(fid,'Model %i (probability %.3f):\n',i,evi(i));
    fprintf(fid,'Estimated value for log10-evidence: %.3f +/- %.3f.\n',results(i).logZ(1)/log(10),results(i).logZ_error/log(10));
    fprintf(fid,misc.titles{1});
    for j=1:length(misc.percentiles_at)
      fprintf(fid,'% .2f    ',misc.percentiles_at(j));
    end
    fprintf(fid,' MaxL@    Mean    +/- dev.\n');
    for j=1:length(models(i).labels)
        fprintf(fid,misc.labels(models(i).labels(j),:));
        for k=1:length(misc.percentiles_at)
          fprintf(fid,ns_print_val(results(i).percentiles(j,k),9));
        end
        fprintf(fid,[ns_print_val(results(i).maxLpar(j),9) ns_print_val(results(i).param_mean(j),9) '+/-' ns_print_val(results(i).param_stddev(j),9) '\n']);
    end
    fprintf(fid,'log10-Maximal likelihood: % .3f\n',results(i).samples(end).logl/log(10));

    if isfield(models(i),'replicate')
      fprintf(fid,'Scaling n:  ');
      for j = 1:length(results(i).prob)
        if isfield(models(i).options,'nlist')
          fprintf(fid,'%-5i ',models(i).options.nlist(j));
        else
          fprintf(fid,'%-5i ',j);
        end
      end 
      fprintf(fid,'\n');
      fprintf(fid,'p-value:    ');
      for j = 1:length(results(i).prob)
        if results(i).prob(j)==0
          fprintf(fid,'0     ');
        elseif results(i).prob(j)==1
          fprintf(fid,'1     ');
        else
          if models(i).options.trackmax<=100
            fprintf(fid,'%.2f  ',results(i).prob(j));
          else
            fprintf(fid,'%.3f ',results(i).prob(j));
          end  
        end
      end 
      fprintf(fid,'\n');
    end
end
fclose(fid);

