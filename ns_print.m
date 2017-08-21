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
      fprintf(fid,ns_print_val(misc.percentiles_at(j),9));
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
      fprintf(fid,'Information content model check:\n');
      fprintf(fid,' Scaling n: ');
      for j = 1:length(results(i).prob)
        if isfield(models(i).options,'nlist')
          fprintf(fid,ns_print_val(models(i).options.nlist(j),8));
        else
          fprintf(fid,ns_print_val(j,8));
        end
      end 
      fprintf(fid,'\n');
      fprintf(fid,' p-value:   ');
      for j = 1:length(results(i).prob)
        fprintf(fid,ns_print_val(results(i).prob(j),8));
      end 
      fprintf(fid,'\n');
      if isfield(models,'checks')
        for j=1:length(models(i).checks)
          fprintf(fid,models(i).checks(j).misc.labels{1});       
          if isfield(models(i).checks(j).misc,'columns')
            fprintf(fid,'\n');
            label2=models(i).checks(j).misc.labels{2};
            fprintf(fid,label2);
            for l=1:length(results(i).checks(j).pvals(1,:));
              fprintf(fid,ns_print_val(models(i).checks(j).misc.columns(l),8));
            end
            fprintf(fid,'\n');
          end
          for k=1:length(results(i).checks(j).pvals(:,1));
            if isfield(models(i).checks(j).misc,'rows')
              fprintf(fid,ns_print_val(models(i).checks(j).misc.rows(k),length(label2)));
            elseif isfield(models(i).checks(j).misc,'columns')
              fprintf(fid,models(i).checks(j).misc.labels{3});
            end
            for l=1:length(results(i).checks(j).pvals(1,:));
              pval=results(i).checks(j).pvals(k,l);
              fprintf(fid,ns_print_val(pval,8));
            end
            fprintf(fid,'\n');
          end
        end
      end
    end
end
fclose(fid);

