function readCalibFiles
global PDR
% function to read calibration data and select appropriate scales and
% atten values for each speaker

pth=PDR.SOUNDS_calib_path;

if(~PDR.DEBUG)
    switch PDR.bird_id
        case 925
            nfiles = 1;
            PDR.SOUNDS_calib_fnames = cell(nfiles,1);
            PDR.SOUNDS_calib_fnames{1} = ['intraural_calib_130322_925D.mat'];
            %PDR.SOUNDS_calib_fnames = {'intraural_calib_130309_925A.mat';...
            %       'intraural_calib_130318_925F.mat'};
            %PDR.SOUNDS_calib_fnames = ['intraural_calib_130222_925B.mat'];
        case 924
            nfiles = 1;
            PDR.SOUNDS_calib_fnames = cell(nfiles,1);
            PDR.SOUNDS_calib_fnames{1} = 'intraural_calib_130309_924B.mat';
            %PDR.SOUNDS_calib_fnames = ['intraural_calib_130221_924A.mat'];
        case 1073
            nfiles = 1;
            PDR.SOUNDS_calib_fnames = cell(nfiles,1);
            PDR.SOUNDS_calib_fnames{1} = ['intraural_calib_130322_1073A.mat'];
            %PDR.SOUNDS_calib_fnames{1} = ['intraural_calib_130309_1073A.mat'];
            %PDR.SOUNDS_calib_fnames = ['intraural_calib_130219_1073A.mat'];
        otherwise
            PDR.SOUNDS_calib_fnames = cell(1,1);
            PDR.SOUNDS_calib_fnames{1}='ERROR';
            warndlg('No calibration data for this bird!!!','!! Warning !!');
            return;
    end
end

tmp = NaN*ones(1,PDR.SOUNDS_num_speakers-1);
PDR.SOUNDS_lead_attens = tmp;
PDR.SOUNDS_speaker_scales_lead = tmp;
PDR.SOUNDS_speaker_scales_lag = tmp;

if(PDR.DEBUG>0)
    tmp=ones(1,PDR.SOUNDS_num_speakers-1);
    PDR.SOUNDS_speaker_scales_lead = 32760.*tmp;
    PDR.SOUNDS_speaker_scales_lag = 32760.*tmp;
    PDR.SOUNDS_lead_attens = 0.*tmp;
    spkr_list=[];
else
    
    expt_locs = [PDR.SOUNDS_elevations; PDR.SOUNDS_azimuths];
    spkr_list = PDR.SOUNDS_speaker_numbers;
    spkr_list = spkr_list(find(spkr_list ~= 1)); % remove lead speaker
    
    % go through each calibration file and grab speaker calibration data
    for j = 1:size(PDR.SOUNDS_calib_fnames,1)
        
        % load calibration data
        load([pth PDR.SOUNDS_calib_fnames{j}]);
        
        calib_locs = CALIB_PDR.locations(:,2:end); % get all locations except lead
        
        % get speaker attenuations & scales
        spkr_IDs = CALIB_PDR.equalized_lead_attens{1,2};
        attens = CALIB_PDR.equalized_lead_attens{2,2};
        speaker_data = CALIB_PDR.speaker_data;
        
        for i0 = spkr_list
            % check if calibration file has calibration data for this speaker
            idx=find(spkr_IDs == i0);
            if length(idx)
                % check if calib data is for same location
                el = expt_locs(1,i0); az = expt_locs(2,i0);
                if el == calib_locs(1,idx) && az == calib_locs(2,idx)
                    % store attenuation for this speaker
                    PDR.SOUNDS_lead_attens(i0-1) = attens(idx);
                    % store scales for this speaker
                    lead_coeffs=speaker_data(idx).coeffs_both_fit_lead_scales;
                    m_lead=lead_coeffs(1,1);
                    b_lead=lead_coeffs(2,1);
                    PDR.SOUNDS_speaker_scales_lead(i0-1) = round(10^((PDR.SOUNDS_SPL-b_lead)/m_lead));
                    lag_coeffs=speaker_data(idx).coeffs_both_fit_lag_scales;
                    m_lag=lag_coeffs(1,1);
                    b_lag=lag_coeffs(2,1);
                    PDR.SOUNDS_speaker_scales_lag(i0-1) = round(10^((PDR.SOUNDS_SPL-b_lag)/m_lag));
                    spkr_list = spkr_list(find(spkr_list ~= i0)); % done with this speaker!
                end
            end
        end
        
    end
end % debug conditions
if length(spkr_list)
    warndlg(['Could not find calib data for speaker numbers: ' num2str(spkr_list)],'!! Warning !!');
    PDR.SOUNDS_calib_fname='ERROR';
end

% check speaker scales to be in range for TDT
scales=[PDR.SOUNDS_speaker_scales_lead PDR.SOUNDS_speaker_scales_lag];
if ~isempty(find((scales>32760)==1))
    warndlg('Some of the scales are too LARGE for the TDT!!!','!! Warning !!');
    PDR.SOUNDS_calib_fname='ERROR';
end
if ~isempty(find((scales<=0)==1))
    warndlg('Some of the scales are too SMALL for the TDT!!!','!! Warning !!');
    PDR.SOUNDS_calib_fname='ERROR';
end