�WPCM  
     �� 2 �       B       R    Z   g        �    #| x          Courier Courier Bold HP LaserJet 4M Plus (PCL)            HL4MPLPC.PRS �x �
   @� ��  � \��KX@�� 2      X  �       K                      <��� ? x x x ,   :W�x �6X   @� �`7�X@��������������������������������������� ? x x x ,  ���x �    `� �7�X�����������������������������������������nel to a file (type 2 or 3)
%
%         filename: string with filename
%         c����el: vector with data
%         index: number of channel (index > 0)       
_______  �  � ?    �� ��\ ��   �function signal = ��abs2db��(samples)
�  � ? �  ���ABS2DB�� changes the complex elements or elements with absolute magnitude tomagnitude with logarithm scaling in [dB].

___________________
�  � ? � ���ANALYSE ��guides through the different visualisation methods of impulse responsesand transfer functions of the external ears.

_________________________
�  � ?  �function hz = ��bark2hz��(bark)
�  � ? � ���BARK2HZ ��changes frequency value from the nonlinear 'bark' - scale to the linearfrequency scale in 'Hz'.

________________________
�  � ? �
 �function f=��checkfft��(H);
�  � ? � ���CHECKFFT�� Check complex transfer function to determine whether H containssamples in the range 0...pi (returns 1) or in the range 0...2*pi (returns 0)

________________________
�  � ? � �function [y,h] = ��cntaxes��(fig)
�  � ? � ���CNTAXES ��Count axes objects.

________________________
�  � ? � �function [y,h] = ��cntlines��(ax)
�  � ? � ���CNTLINES ��Count line axes objects.

________________________
�  � ?  �function [n,y] = ��cntuictl��(fig,style)
�  � ? � ���CNTUICTL ��Count figure uicontrols objects.

________________________
�  � ? 8 ���CONT_ING ��contour plot of the actual signal matrix, that is resampled with         a lower sampling rate for smoothing the contour lines.

________________________
function�  � ?   ���contplot��(file_format,filetype,filename,first_channel,second_channel,first_phase�,second_phase,direction_matrix,n_directions,n_channels,n_lines,xaxis,scaling)
�  � ? � ���CONTPLOT ��plots frequency signals as a contour plot of each channel's magnitudeand phase.

________________________
function[first_channel,second_channel,first_phase,second_phase,n_lines,begin_of_signal,�end_of_signal,xaxis] =�  � ? (# ���cut_sig��(filetype,first_channel,second_channel,first_phase,second_phase,n_direct�ions,n_channels,n_lines,sampling_rate,begin_of_signal,end_of_signal,xaxis)
�  � ? �$ ���CUT_SIG ��cuts a signal in the selected range if cutting is selected.

________________________
�  � ? ' �function ��delay��(filetype,filename,first_phase,second_phase,...
    direction_matrix,n_directions,n_channels,n_lines,xaxis)
�  � ? �( ���DELAY�� plots the phase delay or group delay of signals with the filetype 5 and6.�   h)        0*0*0*��  Ԍ�________________________
�  � ? �  �function ��delaylog��(filetype,filename,first_phase,second_phase,...
     direction_matrix,n_directions,n_channels,n_lines,xaxis)
�  � ? X ���DELAYLOG ��plots the phase delay or group delay of filetype 5 and 6.

________________________
function�  � ? x ���delayout��(filetype,filename,n_o_subplot,headline,n_o_direction,direction_matrix,�xaxis,n_lines)
�  � ?  ���DELAYOUT�� generates a layout for 'delay plots'

________________________
�  � ? `	 �function dy = ��derive��(y,dx);
�  � ? (
 ���DERIVE ��dy=derive(y,dx) computes dy/dx

________________________
�  � ? � �function s=��dinput��(prompt,default,mode);
�  � ? H ���DINPUT�� Keyboard input with default value. s=DINPUT(prompt,default,mode) 
Invokes the MATLAB function INPUT with parameters prompt and mode. If thereturned value from INPUT is empty DINPUT returns the value default.

________________________
�  � ? 0 �function ��dirlabel��(filetype,n_o_direction,direction_matrix)
�  � ? � ���DIRLABEL�� labels the plots with the direction

________________________
�  � ? P �function [] = ��dispdata��(filename);
�  � ?  ���DISPDATA ��displays the measured impulse responses and transfer-functions.
   PRELIMINARY VERSION !!

________________________
�  � ? 8 �function H = ��earfilt��(selected_frq,number_of_lines,resolution,sampling_rate)
�  � ?   ���EARFILT ��evaluates a filter function of the internal ear, with a criticalbandwidth.

________________________
�  � ?   �function ex = ��energy��(x);
�  � ? � ���ENERGY�� returns the energy of vector x. (the summed square of all elements).

________________________
�  � ? @ �function etch = ��etc��(h)
�  � ?   ���ETC�� returns the energy time curve of the (real) vector h.

________________________
�  � ? `" �function [h] = ��fade��(h, b1, e1, b2, e2);
�  � ? (# ���FADE ��Fading function. 

________________________
function [number_of_figures,figure_list,number_of_subplots] =�  � ? H& ���fg_cntrl��(max_plot_in_figure,n_directions,xaxis)
�  � ? ' ���FG_CNTRL ��evaluates the necessary number of figures for plotting and opens them.Returns # of figures, list of handels to each fig, # subplots in each fig. 

________________________�   h)         0*0*0*��  Ԍ�  � ?    �function [h] = ��findax��(fig)
�  � ? �  ���FINDAX�� returns a vector 'h' with the handles of all axes, that are situated inthe figure with the handle 'fig'.

________________________
�  � ? � �function [fg,fm] = ��freqgrp��(fftlen,samplerate);
�  � ? � ���FREQGRP ��Upper cutoff frequencies fg and center frequencies fm of Zwickersfrequency groups.

________________________
function selected_frq =�  � ? � ���frqcntrl��(frq,n_lines,begin_of_signal,end_of_signal,xaxis)
�  � ? `	 ���FRQCNTRL ��returns the value of the sampled frequency, that is the nearest to thechoosen frequency, if the choosen one is not sampled exactly.

________________________
�  � ? � �function scaling = ��frqscale��()
�  � ? H ���SCALING ��returns a string 'scaling' defining the kind of scaling for thefrequency axis. The scaling can be linear ('lin') or logarithm ('log')

________________________
�  � ? h �function Hs = ��fsmooth��(H,width,fftlen)
�  � ? 0 ���FSMOOTH�� octave ore third octave smoothing of transferfunction H.

________________________
�  � ? � �function result = ��gauss��(value)
�  � ? P ���GAUSS�� returns the next integer value greater than the given number 'value'.

________________________
�  � ? � �function [sample,samplestr] =�� getsample��(mode)
�  � ? p ���GETSAMPLE ��return sampling points for monaural and interaural optimization.

________________________
�  � ? � �function ��grcntrl��()
�  � ? � ���GRCNTRL�� User menu for selecting plot windows.

________________________
�  � ? � �function ��grcntrl_��(id)
�  � ? � ���GRCNTRL_ ��Is called from GRCNTRL.

________________________
�  � ?   �function [hout,t] = ��h_max_en��(h,N)
�  � ? �  ���H_MAX_EN�� A rectangular window of length N is shifted over the impulse responseh and the energy inside the window is computed as a function of the startposition of the window. The point t where the energy the first time becomes amaximum is detected. A modified impulse response hout is generated by ignoringall samples in front of that point.

______________________
�  � ? H& �function bark = ��hz2bark��(hz)
�  � ? ' ���HZ2����BARK�� changes frequency value from linear 'Hz' - scale to the nonlinear'bark'.

�   h)         0*0*0*��  �________________________
�  � ? �  ���IA_CNTRL�� provides a menu for different kind of plots of interaural level
   differences versus interaural time differences.

________________________
�  � ? � �function ��ia_dir��(ia_level,ia_time,direction_matrix,n_directions)
�  � ? � ���IA_DIR�� labels the directions in a 2D plot of the interaural level differenceversus the interaural time difference.

________________________
function [level_difference,time_difference,frq_axis] = ...
�  � ? � ���ia_eval1��(first_channel,second_channel,n_directions,n_lines,sampling_rate,begin_�of_signal,end_of_signal)
�  � ? (
 ���IA_EVAL1�� evaluates the interaural level and time difference for all givendirections and a single frequency.

________________________
function [level_difference,time_difference,frq_axis] =�  � ?  ���ia_evalx��(first_channel,second_channel,n_directions,n_lines,sampling_rate,begin_�of_signal,end_of_signal,xaxis)
�  � ? � ���IA_EVALX�� evaluates the ILD and ITD for all directions and all frequencies in aselected range. The frequency range must be selected in [bark].

________________________
function [level_grid,level_color,level_contour,azi_interp,...
�  � ? � � elev_interp] = ��ia_fit��(difference,direction_matrix,n_directions,type)
�  � ? P ���IA_FIT ��linear interpolation of the ILDs or ITDs of the measured directions todirections, that are not measured, to get a better resolution.

________________________
�  � ? p �function ��ia_globe��(difference,direction_matrix,frq_axis,type)
�  � ? 8 ���IA_GLOBE�� plots a sphere with a combined plot of the iso - ILDs or ITDs for asingle frequency

________________________
�  � ? X ���IA_LYSE�� guides through the different visualisation methods of ILDs and ITDs.

________________________
�  � ? � �function ��ia_plot��(level_difference,time_difference,direction_matrix,...
      frq_axis,single_plane)
�  � ? @ ���IA_PLOT�� plots for a horizontal plane the ITD and ILD about the changing azimuthangle for a single frequency.

________________________
function�  � ? (# ���ia_plot1��(level_difference,time_difference,direction_matrix,n_directions,frq_axi�s,single_plane)
�  � ? �$ ���IA_PLOT1�� plots the ILD vs ITD

�   H&         0*0*0*��  �________________________
function�  � ? � ���ia_plot2��(level_difference,time_difference,direction_matrix,n_directions,frq_axi�s)
�  � ?   ���IA_PLOT2�� plots a 3D plot of the ILD versus ITD for all evaluated frequencies.The different directions are plotted in different colors.

________________________
�  � ? @ �function [selected_frq,selected_line] = ��ia_s_frq��(frq_axis,n_iaural)
�  � ?  ���IA_S_FRQ�� selects a frequency for plotting interaural signal differences.

________________________
function�  � ? (
 ���ia_save��(save_name,level_difference,time_difference,direction_matrix,n_direction�s,sampling_rate,frq_axis,comment)
�  � ? � ���IA_SAVE�� saves the evaluated ILDs and ITDs in the MTLSIG - format.

________________________
�  � ?  �function ��irplot��(h,samplerate,gtitle,color,mintime,maxtime);
�  � ? � ���IRPLOT�� Plot one or more real (time) signals.

________________________
function [level_grid,level_color,level_contour,azi_interp,elev_interp] =�  � ? � ���iso_fit��(magn,direction_matrix,n_directions,min_level,max_level,level)
�  � ? � ���ISO_FIT�� linear interpolation of the magnitudes of the measured directions todirections, that are not measured, to get a better resolution.

________________________
function�  � ? � ���iso_mesh��(isofrequency,level_grid,level_color,level_contour,azi_interp,elev_inte�rp)
�  � ? 8 ���ISO_MESH�� plots lines of the same level in a 3D-meshplot. The different levelsare displayed with a pseudocolor plot and additional contour lines.

________________________
function�  � ?   ���isoglobe��(isofrequency,level_grid,level_color,level_contour,azi_interp,elev_inte�rp)
�  � ? � ���ISOGLOBE�� plots lines of the same level on a sphere. Different levels aredisplayed with a pseudocolor plot and additional contour lines.

________________________
�  � ? �  ���ISOLYSE�� guides through the visualisation methods for isolines of transferfunctions of the external ears.

________________________
�  � ? �# �function ��isoplot��(filename,first_channel,second_channel,selected_channel,...direction_matrix,n_directions,n_lines,sampling_rate,begin_of_signal,end_of_sign�al,plot_kind)
�  � ? H& ���ISOPLOT�� plots a 3D plot of isolevels at a selected frequency. Two different

�   �'         0*0*0*��  �________________________
�  � ? �  �function okay=��isstable��(p);
�  � ? � ���ISSTABLE��  Check polynom with coefficients in p for roots outside or on the unitcircle. 

________________________
�  � ? � �function ��l_contr��(cs)
�  � ? x ���L_CONTR�� adds height labels to the current contour plot using the contourstructure 'cs', the returning output from the CONTOUR routine. 

________________________
�  � ? � �function [K,ec] = ��lattice��(A)
�  � ? `	 ���LATTICE�� [K,ec]=LATTICE(A) calculates lattice reflection coefficients K frompolynomial coefficents given in A.

________________________
�  � ? � �function level_difference = ��leveldif��(left_X,right_X)
�  � ? H ���LEVELDIF ��evaluates the interaural level difference by the following method.  The transfer functions of the external ears are filtered by a ROEX-filter (seeEARFILT) with the critical bandwidth depending on the choosen frequency tofilter. The intensities of the the left and the right channel are evaluated andthe level difference related to the left ear is evaluated and returned in [dB].

________________________
�  � ? � �function ��lsequal2��(displaymode,listname)
�  � ? � ���LSEQUAL2 ��Compute random phase signals with loadspeaker equalization for theloudspeakers for HRTF measurement.

________________________
�  � ? � �function hi=��lsqinv��(h,N,L,wlp)
�  � ? p ���LSQINV�� Least square inverse filter design

________________________
function�  � ? � ���m_layout��(file_format,filetype,filename,n_o_subplot,headline,n_o_direction,direc�tion_matrix,xaxis,n_lines)
�  � ?   ���M_LAYOUT ��generates a layout for plots of the magnitude of transfer functions.

________________________
function�  � ? @ ���magn_log��(file_format,filetype,filename,first_channel,second_channel,direction_m�atrix,n_directions,n_channels,n_lines,xaxis)
�  � ? �  ���MAGN_LOG�� Plots the signals magnitude (first_channel,second_channel) over'xaxis' with a logarithm x-axis.

________________________
�  � ? �# �function ��meantf��(filelist,outname,t)
�  � ? �$ ���MEANTF ��computes the mean magnitude transfer function for the specified filelist and writes it to the file outname.

________________________
�  � ? �' ���MESH_LAY ��generates a layout for 'mesh plots'. Different layouts are providedfor the different file formats.
�   h)         0*0*0*��  Ԍ�              0*0*0*��  �________________________
�  � ? �  ���MESHING ��Meshplot of the actual signal matrix. Sliders for changing the point ofview are provided, too.

________________________
function�  � ? � ���meshplot��(file_format,filetype,filename,first_channel,second_channel,...first_phase,second_phase,direction_matrix,n_directions,n_channels,n_lines,xaxis�,scaling)
�  � ?  ���MESHPLOT�� generates a meshplot of the magnitude for each channel. For frequencysignals the phase for each channel is plotted, too.

________________________
�  � ? (
 �function out = ��mirror��(n_lines,in)
�  � ? �
 ���MIRROR�� conjugate complex mirroring of the matrix 'in'. 'n_lines' is the numberof columns of the matrix. If n_lines is even the resulting output matrix has 2* n_lines columns, if n_lines is odd the resulting output matrix, has 2 *(n_lines-1) columns.

________________________
�  � ? � �function [H1, h1, H2, ofs] = ��mklsdata��(N, M, L, h)
�  � ? h ���MKLSDATA�� make data arrays for least-square approximation

________________________
�  � ? � �function [H, h] = ��mklsrdat��(N, M, L, r, zlag)
�  � ? � ���MKLSRDATA ��make data arrays for least-square approximation with autocorrelationdata

________________________
�  � ? � �function sample=��mksample��(fg,N,samplerate,fftlen);
�  � ? p ���MKSAMPLE  ��Return nonuniformly spaced sampling points.

________________________
�  � ? � �function ��mksinsi1��(filename,magnitude,frequency);
�  � ? � ���MKSINSI1  ��writes a sine signal to the specified file using SI1 format. 

________________________
�  � ? � �function sample = ��mmsample��(H,minsample,maxsample);
�  � ? � ���MMSAMPLE ��Compute sampling points near the maxima and minima points oflog10(abs(H)). sample does contains no values below MINSAMPLE or aboveMAXSAMPLE.

________________________
�  � ? �! ���MO_CNTRL�� is a submenu called from 'analyse', for the control of motions. Motion always displays signal's magnitude. 

________________________
�  � ? �$ �function fps = ��mo_speed��()
�  � ? �% ���MO_SPEED ��changes the speed of motions. 

________________________
function mtlcont(data,azi,ele,n)
plots contour lines on sphere of mtlplots
�   h)         0*0*0*��  Ԍ�     	         0*0*0*��  �________________________
function mtlcorhe(filename)
corrects header information

________________________
�  � ? � �function ��mtlcut��(infile, outfile, lwindow, n_lines_out)
�  � ? � ���MTLCUT�� cuts the most important part of an impulse response.

________________________
�  � ?  �function ��mtldiv��(file1, file2, outfile)


________________________
�  � ? (
 �function ��mtlfade��(infile,outfile)
ramps on and off at 3 & 10 kHz

________________________
�  � ? H �function index = ��mtlfd��(dir,azimuth, elevation)


________________________
�  � ? h �function ��mtlfft��(infile,outfile,len)
converts time file into frequency file

________________________
�  � ? � �function ��mtlia��(file1, outfile)


________________________
�  � ? � �function ��mtlifft��(infile,outfile)
changes frequency file into time file

________________________
�  � ? � �function ��mtlinv��(infile,outfile,f1, f2 ,f3 ,f4)


________________________
�  � ? � �function ��mtllow��(infile, outfile, frequency)


________________________
�  � ?   �function ��mtlmean��(infile,outfile)

________________________
�  � ? `" �function ��mtlmerge��(file_list, outputfile)
merges all files in list into a cataloque

________________________
�  � ? �% �function ��mtlplot��(catalogue, elevation, azimuth)
�  � ? H& ���MTLPLOT�� plots the impulse response or transfer function of the selecteddirection or channel

�   �(
         0*0*0*��  �________________________
function mtlproc(listname)


________________________
�  � ? � �function channel=��mtlrch��(filename, index);
�  � ? � ���MTLRCH�� reads one channel from a file (type 2 or 3)

________________________
�  � ?  �function delay=��mtlrdel��(filename,n_channels)
�  � ? � ���MTLRDEL ��reads delay from MTLSIG catalogue

________________________
�  � ? (
 �function direction_matrix=��mtlrdir��(filename,n_directions)
�  � ? �
 ���MTLRDIR ��reads direction matrix of a mtlsig file

________________________
�  � ? H �function direction_matrix=��mtlrdir��(filename,n_directions)


________________________
function [filetype,info_blocks,n_channels,n_lines,sampling_rate,first_line,last_line,n_directions,comment1, comment2] =�  � ? � ���mtlrh��(filename)
�  � ? � ���MTLRH ��reads the header of a SIGMAP and MTLSIG file.

________________________
�  � ?  �function ��mtlscale��(infile, outfile, factor)


________________________
�  � ? 8 �function status=��mtlseek��(fid, offset, origin)


________________________
�  � ? X �function ��mtlsf��(infile, outfile, new_sampling_rate)


________________________
�  � ? x �function filetype = ��mtlfiletype��(filename)


________________________
�  � ? �! �function status = ��mtlwch��(filename,channel, index)
MTLWCH writes one channel to a file (type 2 or 3)

________________________
�  � ? �$ �function status = ��mtlwdel��(filename,delay_l,delay_r)

________________________
�  � ? ' �function status = ��mtlwdel2��(filename,delay)
MTLWDEL writes delay to catalogue

�   h)         0*0*0*��  �________________________
�  � ? �  �function status=��mtlwdir��(filename,direction_matrix)
MTLWDIR writes direction matrix of a mtlsig file

________________________
�  � ? � �function  ��mtlwh��(filename,filetype,info_blocks,n_channels,n_lines,sampling_rate,...
first_line,last_line,n_directions,comment1, comment2)
�  � ? @ ���MTLWH ��writes the header of a SIGMAP and MTLSIG file.

________________________
functionp_layout(file_format,filetype,filename,n_o_subplot,headline,n_o_direction,direc�tion_matrix,xaxis,n_lines)
�  � ? �
 ���P_LAYOUT�� generates a layout for the plots of the phase of transfer functions.

________________________
function�  � ?  ���ph_log��(file_format,filetype,filename,first_channel,second_channel,first_phase,s�econd_phase,direction_matrix,n_directions,n_channels,n_lines,xaxis)
�  � ? � ���PH_LOG ��Plots the phase, phase delay or group delay of the signals over 'xaxis'with a logarithm x-axis scaling. 

________________________
�  � ? � �function ��ph_plot��(file_format,filetype,filename,first_channel,second_channel,...
first_phase,second_phase,direction_matrix,n_directions,n_channels,n_lines,xaxis�)
�  � ?  ���PH_PLOT�� Plots the phase, phase delay or group delay of the signals over 'xaxis'with a linear x-axis scaling.

________________________
�  � ? 8 �function ��phasplot��(H,samplerate,gtitle,color,minfreq,maxfreq,correct,unwraptol);
�  � ?   ���PHASPLOT�� Plot phase of one or more complex signals in degrees.

________________________
�  � ? X �function single_plane = ��planedct��(direction_matrix)
�  � ?   ���PLANEDCT�� detects if the direction matrix includes a single horizontal plane.Then 'single_plane' includes the string 'elev', else the string 'none' isreturned.

________________________
�  � ?   �function [so,nc]=��pltstyle��(si);
�  � ? �  ���PLTSTYLE�� Make correct plot style string

________________________
�  � ? (# ���POLARNET�� plots a net of lines, like a globe on an already existing sphere, forbetter spatial orientation. POLARNET is based on the function 'sphere', butuses a higher resolution, for better drawing approximation.

�   H&         0*0*0*��  �________________________
function�  � ? � ���psc_con��t(file_format,filetype,filename,first_channel,second_channel,first_phase�,second_phase,direction_matrix,n_directions,n_channels,n_lines,xaxis,scaling)
�  � ?   ���PSC_CONT ��plots for frequency signals a combined pseudocolor and contour plot ofeach channels magnitude and phase.

________________________
�  � ? @ ���PSC_ING�� Pseudocolor and contour plot of the actual signal matrix, that isresampled with a lower sampling rate for smoothing the contour lines.

________________________
function[file_format,first_channel,second_channel,first_phase,second_phase,direction_ma�trix,n_directions,n_channels,n_lines,sampling_rate,begin_of_signal,end_of_signa��  � ? � �l,xaxis,comment] = ��r_frq��(fid,filetype)
�  � ? � ���R_FRQ�� Reads the specific frequency signal attributes. 

________________________
function[file_format,first_channel,second_channel,direction_matrix,time_delay,...n_directions,n_channels,n_lines,sampling_rate,begin_of_signal,end_of_signal,xax��  � ? 0 �is,comment] = ��r_time��(fid,filetype)
�  � ? � ���R_TIME ��Reads the specific time signal attributes.

________________________
function[file_format,filetype,first_channel,second_channel,first_phase,second_phase,dir�ection_matrix,time_delay,n_directions,n_channels,n_lines,sampling_rate,begin_of��  � ? � �_signal,end_of_signal,xaxis,comment] = ��read��(filename)
�  � ? p ���READ�� Reads catalog of hrtf or hrir signals

________________________
�  � ? � �function direction_matrix=��read_dir��(filename,n_directions)


________________________
function[isoloaded,criterium,filetype,sampling_line,value,delta,iso_points,iso_array] =�  � ? x ���read_iso��(filename,isoloaded,criterium,filetype,sampling_line,value,delta,iso_po�ints,iso_array)
�  � ?   ���READ_ISO ��reads the data of an iso-file

________________________
�  � ? `" �function [h,attr]=��readdata��(catname,azimuth,elevation,channel,mode,display);
�  � ? (# ���READDATA��  Read measured impulse response from catalogue.

________________________
�  � ? �% �function [d,mindelay]=��readdel��(catname)
�  � ? H& ���READDEL��  returns the delay times and the minimum delay time of the impulseresponses in the specified catalogue.

�   �(         0*0*0*��  �________________________
�  � ? �  �function [d,ndirs]=��readdirs��(catname)
�  � ? � ���READDIRS�� returns the available directions and the number of directions storedin the specified catalogue.

________________________
function [filetype,info_blocks,n_channels,n_lines, sampling_rate,first_line,�  � ? x �last_line,n_directions,comment1, comment2] = ��readh��(filename)
�  � ? @ ���READ_H�� Reads the header of a SIGMAP and MTLSIG file.

________________________
�  � ? � �function [h,attr]=��readhp��(catname,channel,mode,display);
�  � ? `	 ���READHP��  Read measured headphone impulse response for given catalogue.

________________________
�  � ? � �function list=��readlist��(cmdfile)
�  � ? � ���READLIST�� reads a ascii file containing a file list

________________________
�  � ? � �function [h,attr]=��readls��(catname,number,mode,display);
�  � ? � ���READLS��  Read measured loudspeaker impulse response from catalogue.

________________________
�  � ? � �function [Hl,Hr,azi,ele]=��readsc2��(filename);
�  � ? � ���READSC2�� reads SC2 file

________________________
�  � ?  �function [h,header]=��readsi1��(filename);
�  � ? � ���READSR1�� reads SI1 file

________________________
�  � ? 8 �function [hl,hr]=��readsig2��(filename);
�  � ?   ���READSIG2�� reads SIG2 file

________________________
�  � ? X �function [h,header]=��readsr1��(filename);
�  � ?   ���READSR1�� reads SR1 file

________________________
�  � ? x �function [hl,hr]=��readsr2��(filename);
�  � ? @ ���READSR2�� reads SR2 file

________________________

�  � ? `" �function X = ��rfft��(x,len)


________________________
�  � ? �% �function fa=��samplfrq��()
SAMPLEFRQ returns the default sampling frequency.
used from nearly all other functions, so the samplerate could be changedeasily.

�   h)         0*0*0*��  �________________________
�  � ? �  �function ch=��seldirs��(dirs,ch);
SELDIRS displays the available directions (44!) and enables the user to selectspecial directions.

________________________
�  � ? � �function y=��smooth��(x,M);
�  � ? x ���SMOOTH�� Smooth data vector x using M values on the left and right of each point.

________________________
�  � ? � �function ��sph_cont��(level_contour)
�  � ? � ���SPH_CONT ��plots contour lines on a sphere and labels them. 'level_contour' isthe matrix C returned by the function 'contourc'.

________________________
�  � ? � �function a=��sphdist��(sp1,sp2)
�  � ? � ���SPHDIST ��calculates the distance of two points on a sphere

________________________
�  � ? � �function ��sqrplot��()
�  � ? � ���SQRPLOT�� produces equal sized axis for a graph.

________________________
function�  � ? � ���t_layout��(file_format,filetype,filename,n_o_subplot,headline,n_o_direction,direc�tion_matrix,xaxis,n_lines)
�  � ? P ���T_LAYOUT�� generates a layout for plots of impulse responses. Different layoutsare available, depending on the different file formats.

________________________
�  � ? p �function grpdelay = ��tgroup��(H,samplerate);
�  � ? 8 ���TGROUP�� Group delay for vector H.

________________________
�  � ? � �function ��tgrpplot��(H,samplerate,gtitle,color,minfreq,maxfreq);
�  � ? X ���TGRPPLOT ��Plot group delay of one or more complex signals.

________________________
�  � ? � �function time_difference = ��timediff��(frq,n_lines,sampling_rate,left_X,right_X)
�  � ? x ���TIMEDIFF ��evaluates the ITD by the following method.
   The transfer functions of the external ears are filtered by a ROEX-filter(see EARFILT) with the critical bandwidth depending on the choosen frequency tofilter.
For frequencies less than 1.6 kHz the displacement of the carrier frequenciesis evaluated. This is done by multiplication of the left external ear'stransfer function with the conjugate complex transfer function of the rightear. The resulting single band cross power density is conjugate complexmirrored. The inverse.  Fourier transformed leads to the cross correlation ofthe two
impulse responses. The displacement of the maximum from the middle describesthe interaural time difference. 

�   �(         0*0*0*��  �________________________
�  � ? �  �function ��tphplot��(H,samplerate,gtitle,color,minfreq,maxfreq,correct,unwraptol);
PRELIMINARY
�  � ? X ���TPHPLOT�� Plot phase delay of one or more complex signals in degrees.

________________________
function [u_s_signal,new_xaxis,newlines] =�  � ? x ���u_sample��(signal,xaxis,n_directions,n_lines,resample_factor)
�  � ? @ ���U_SAMPLE�� resamples the in the matrix 'signal' included data with a lower
sampling rate after lowpass filtering.

________________________
�  � ? `	 �function gridhandle =�� usergrid��(X,Y,C);
�  � ? (
 ���USERGRID  ��Draw grid lines 

________________________
�  � ? � �function e=��winenerg��(h,wl,n1,n2)
�  � ? H ���WINENERG�� Normalized energy in rectangular window e=winenerg(h,wl,n1,n2)

________________________
�  � ? � �function  ��writeh��(filename,filetype,info_blocks,n_channels,n_lines,sampling_rate,first_line,last_line,n_directions,comment1, comment2)
�  � ? 0 ���WRITEH ��writes the header of a SIGMAP and MTLSIG file.

________________________
�  � ? � �function []=��writesi1��(data,filename,header);
�  � ? P ���WRITESI1  ��Write data to file using SI1 format.

________________________
�  � ? � �function ��writesr2��(filename,hl,hr);
�  � ? p ���WRITESR2�� writes SR2 file with given filename, data hl for left channel and datahr for right channel.

________________________
�  � ? � �function ��writetxt��(fname,str)
�  � ? X ���WRITETXT ��writes a text string to a file.

________________________
�  � ? � �function [ytic,lab,ang1] = ��y_axis��(direction_matrix,n_directions)
�  � ? x ���Y_AXIS�� initializes a vector 'ytic' with the angles, that have the greatestrange. 'lab' includes a string matrix with the angles for labeling. 'ang1'provides which angle has the greatest range.

________________________