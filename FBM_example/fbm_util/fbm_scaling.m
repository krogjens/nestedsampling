function scaled_steps = fbm_scaling(steps,n)

if n> 1
   transpose = 0;
   if size(steps,2) < size(steps,1)
      steps = steps';
      transpose = 1;
   end
   for i = 1:floor(length(steps)/n)
      indexes = n*(i-1) + (1:n);
      scaled_steps(:,i) = sum(steps(:,indexes)')';
   end
   if transpose == 1;
      scaled_steps = scaled_steps';
   end
else
   scaled_steps = steps;
end
