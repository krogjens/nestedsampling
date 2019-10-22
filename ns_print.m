function ns_print(results,models,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prints out a summary of 'results' to the text file ['misc.data_id',misc.nssummary]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = [misc.data_id,misc.nssummary];
if isfield(misc,'append')
  fid = fopen(fname,'a');
  fprintf(fid,misc.append);
else
  fid = fopen(fname,'w');
end

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
    perc_text='Percentile:';
    for j=1:min(length(misc.labels(1,:)),length(perc_text))
      fprintf(fid,perc_text(j));
    end
    for j=length(perc_text):(length(misc.labels(1,:))-1)
      fprintf(fid,' ');
    end
    for j=1:length(misc.percentiles_at)
      fprintf(fid,ns_print_val(misc.percentiles_at(j),9));
    end
    fprintf(fid,' MaxL@    Mean    +/- dev.\n');
    for j=1:min(length(models(i).labels),length(results(i).param_mean))
      if models(i).labels(j)>0
        fprintf(fid,misc.labels(models(i).labels(j),:));
        for k=1:length(misc.percentiles_at)
          fprintf(fid,ns_print_val(results(i).percentiles(j,k),9));
        end
        fprintf(fid,[ns_print_val(results(i).maxLpar(j),9) ns_print_val(results(i).param_mean(j),9) '+/-' ns_print_val(results(i).param_stddev(j),9) '\n']);
      end
    end
    fprintf(fid,'log10-Maximal likelihood: % .3f\n',results(i).samples(end).logl/log(10));

    if isfield(models(i),'replicate')
      fprintf(fid,'Information content model check:\n');
      fprintf(fid,' Scaling n:');
      for j = 1:length(results(i).prob)
        if isfield(models(i).options,'nlist')
          fprintf(fid,ns_print_val(models(i).options.nlist(j),8));
        else
          fprintf(fid,ns_print_val(j,8));
        end
      end 
      fprintf(fid,'\n');
      fprintf(fid,' p-value:  ');
      for j = 1:length(results(i).prob)
        fprintf(fid,ns_print_val(results(i).prob(j),8));
      end 
      fprintf(fid,'\n');
% Print the results of any user defined model checks
      if isfield(models,'checks')
        for j=1:length(models(i).checks)
          fprintf(fid,models(i).checks(j).misc.labels{1});       
          if isfield(models(i).checks(j).misc,'columns')
            fprintf(fid,'\n');
            if length(models(i).checks(j).misc.labels)>1
              label2=models(i).checks(j).misc.labels{2};
              len1stC=length(sprintf(label2));
              fprintf(fid,label2);
            else
              if isfield(models(i).checks(j).misc,'rows')
                fprintf(fid,' R  \\  C');
                len1stC=8;
              else
                fprintf(fid,' Input:  ');
                len1stC=9;
              end
            end
            for l=1:length(results(i).checks(j).pvals(1,:));
              fprintf(fid,ns_print_val(models(i).checks(j).misc.columns(l),8));
            end
            fprintf(fid,'\n');
          end
          for k=1:length(results(i).checks(j).pvals(:,1));
            if isfield(models(i).checks(j).misc,'rows')
              fprintf(fid,ns_print_val(models(i).checks(j).misc.rows(k),len1stC));
            elseif isfield(models(i).checks(j).misc,'columns')
              fprintf(fid,' p-value:');
              for l=1:(len1stC-9);
                fprintf(fid,' ');
              end
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
fprintf(fid,'\n');
fclose(fid);
end

%%% A function that prints a value with a format depending on its magnitude %%%
function txt = ns_print_val(val,len)
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
end
