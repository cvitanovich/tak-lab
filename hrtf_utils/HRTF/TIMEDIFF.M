function time_difference = timediff(frq,n_lines,sampling_rate,...
                                    left_X,right_X)

% function time_difference = timediff(frq,n_lines,sampling_rate,...
%                                     left_X,right_X)
%
% TIMEDIFF evaluates the interaural time difference by the following method.
%          The transfer functions of the external ears are filtered by a 
%          ROEX-filter (see EARFILT) with the critical bandwidth depending on
%          the choosen frequency to filter.
%          For frequencies less than 1.6 kHz the displacement of the carrier
%          frequencies is evaluated. This is done by multiplication of the
%          left external ear's transfer function with the conjugate complex
%          transfer function of the right ear. The resulting single band
%          cross power density is conjugate complex mirrored. The inverse
%          Fourier transformed leads to the cross correlation of the two
%          impulse responses. The displacement of the maximum from the middle
%          describes the interaural time difference. 
%          For frequencies greater than 1.6 kHz the displacement of the
%          envelopes is evaluated. This is done by multiplication of the left
%          external ear's transfer function with the conjugate complex
%          transfer function of the right ear. To get the analytic signal
%          zeros are added. The absolute magnitude of the inverse Fourier
%          transformed leads to the cross correlation of the envelopes. The 
%          displacement of the maximum from the middle describes the
%          interaural time difference. 
%
%          Description of the input arguments:
%          - frq: frequency for which the interaural time difference should be
%                 evaluated (see separation above)
%          - n_lines: number of lines of the filtered transfer functions;
%                     used for mirroring and generating the analytic signal
%          - sampling_rate
%          - left_X: the filtered transfer function of the left ear
%          - right_X: the filtered transfer function of the right ear
%
%          For further information see: leveldif, ia_eval1, earfilt
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


spec = left_X .* conj(right_X);

if frq < 1600
  spec = mirror(n_lines,spec);
  ccf = real(ifft(spec));
else
  spec = [4*spec,zeros(1,n_lines-2)];
  cross_correlation = ifft(spec);
  ccf = abs(cross_correlation);
end;

%figure;
%plot(abs(spec));
%subplot(2,2,4);
%plot(unwrap(angle(spec)));



% evaluate time_difference

zero_pos = (length(ccf))/2;
shift_ccf = [];
shift_ccf(1:zero_pos) = ccf(zero_pos+1:length(ccf));
shift_ccf(zero_pos+1:length(ccf)) = ccf(1:zero_pos);

%figure;
%plot(shift_ccf);
%pause;

[maximum,maxpos] = max(real(shift_ccf));
time_pos = zero_pos - maxpos;
time_difference = time_pos / sampling_rate;

%disp(['maxpos: ',num2str(maxpos)]);
%disp(['timedifference: ',num2str(time_difference)]);