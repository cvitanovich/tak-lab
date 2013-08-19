function time_difference = timediff(frq,n_lines,sampling_rate,...
                                    left_X,right_X)


if frq < 1600
  left_X = mirror(n_lines,left_X,'conj');
  right_X = mirror(n_lines,right_X);
  spec = left_X .* right_X;
  cross_correlation = ifft(spec);
  [maximum,maxpos] = max(cross_correlation);
  time_difference = (maxpos-1) / sampling_rate;
else
  spec = left_X .* right_X;
  % construct analytic signal
  spec = [4*spec,zeros(1,n_lines-2)];
  cross_correlation = ifft(spec);
  abs_ccf = abs(cross_correlation);
  [maximum,maxpos] = max(abs_ccf);
  time_difference = (maxpos-1) / sampling_rate;
end;