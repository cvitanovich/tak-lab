function [result_mean,result_std] = ...
   proc_test_space(datamatrix,...
   location_matrix,...
   time_window,...
   reps_to_include)

disp('Starting data processing...')
rowindex = find(datamatrix(:,1) >= time_window(1) & datamatrix(:,1) <= time_window(2));
data_restrict = datamatrix(rowindex,:);

%Count spikes from appropriate repetitions
disp('Counting spikes for specified repetitions...')
if(isempty(reps_to_include))
   reps_to_include = 1:max(data_restrict(:,2));
end
rowindex = [];
for numinclrep = 1:length(reps_to_include)
   rowindex = find(data_restrict(:,2) == reps_to_include(numinclrep));
   for numcond = 1:length(location_matrix)
      numspikes(numinclrep,numcond) = ...
         length(find(data_restrict(rowindex,3) == location_matrix(1,numcond) & ...
         data_restrict(rowindex,4) == location_matrix(2,numcond)));
   end
   disp(['Total # spikes for rep ' num2str(reps_to_include(numinclrep)) ' = ' ...
         num2str(sum(numspikes(numinclrep,:,:)))])
end

%Calculate mean and standard deviation
disp('Calculating mean and standard deviation...')
if(size(data_restrict,2) < 4)
   result_mean = mean(numspikes,1);
   result_std  = std(numspikes,0,1);
elseif(size(data_restrict,2) == 4)
   result_mean = squeeze(mean(numspikes,1));
   result_std  = squeeze(std(numspikes,0,1));
end

return