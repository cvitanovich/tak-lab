%IAMultRegress: GUI to control a Multiple Regression Session for ILDAlone data


%*******************************************************************************
%	The Graphical User Interface
%*******************************************************************************

pos_incr = 20;
%Figure window
iamrfig_wd = 300; iamrfig_ht = 500;
hiamrfig = figure('Units','pixels',...
   'Position',[300 0 iamrfig_wd iamrfig_ht],...
   'Name','DATATOOL: ILDAlone MultRegress',...
   'NumberTitle','off');
tot_ht = iamrfig_ht;
ind_ht = 20; tot_ht = tot_ht - ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 tot_ht iamrfig_wd ind_ht],...
   'ForegroundColor','blue',...
   'FontWeight','bold',...
   'String', 'DATATOOL: ILDAlone Multiple Regression');

%Specify information for multiple regression
%Datafile
ind_ht = 20; tot_ht = tot_ht - ind_ht - 10;
hdatafilehdr = uicontrol('Parent',hiamrfig,...
   'Style','text',...
   'Position',[0 tot_ht iamrfig_wd ind_ht],...
   'String', 'DataFile:');
ind_ht = 20; tot_ht = tot_ht - ind_ht;
hdatafile = uicontrol('Parent',hiamrfig,...
   'Style','text',...
   'Position',[0 tot_ht iamrfig_wd ind_ht],...
   'String', [XStimParams.datadir XStimParams.datafile]);
ind_ht = 20; tot_ht = tot_ht - ind_ht;
uicontrol('Parent',hiamrfig,...
   'Style','pushbutton',...
   'Units','pixels',...
   'Position',[0 tot_ht 100 ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select DataFile',...
   'Callback','datafilediagbox');

%HRTFfile
ind_ht = 20; tot_ht = tot_ht - ind_ht - 10;
hhrtffilehdr = uicontrol('Parent',hiamrfig,...
   'Style','text',...
   'Position',[0 tot_ht iamrfig_wd ind_ht],...
   'String', 'HrtfFile:');
ind_ht = 20; tot_ht = tot_ht - ind_ht;
hhrtffile = uicontrol('Parent',hiamrfig,...
   'Style','text',...
   'Position',[0 tot_ht iamrfig_wd ind_ht],...
   'String', XStimParams.ildalonefiltfile);
ind_ht = 20; tot_ht = tot_ht - ind_ht;
uicontrol('Parent',hiamrfig,...
   'Style','pushbutton',...
   'Units','pixels',...
   'Position',[0 tot_ht 100 ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select HrtfFile',...
   'Callback','hrtffilediagbox');


%Testnumber for ILDAlone Measurement
ind_ht = 20; tot_ht = tot_ht - ind_ht - 10;
uicontrol('Parent',hiamrfig,...
   'Style','text',...
   'Position',[0 tot_ht 200 ind_ht],...
   'String', 'Test Number for ILDAlone Measurement');
hiatestnum = uicontrol('Parent',hiamrfig,...
   'Style','edit',...
   'Position',[210 tot_ht 50 ind_ht],...
   'String', num2str(XStimParams.testnum));

%Testnumber for ILDFreq Measurement
ind_ht = 20; tot_ht = tot_ht - ind_ht;
uicontrol('Parent',hiamrfig,...
   'Style','text',...
   'Position',[0 tot_ht 200 ind_ht],...
   'String', 'Test Number for ILDFreq Measurement');
hiftestnum = uicontrol('Parent',hiamrfig,...
   'Style','edit',...
   'Position',[210 tot_ht 50 ind_ht],...
   'String', num2str(XStimParams.testnum));

%Set independent variable for regression
ind_ht = 20; tot_ht = tot_ht - ind_ht - 30;
uicontrol('Parent',hiamrfig,...
   'Style','text',...
   'Position',[0 tot_ht 200 ind_ht],...
   'String', 'Regress on ILDFreq Data');
hregressif = uicontrol('Parent',hiamrfig,...
   'Style','checkbox',...
   'Position',[210 tot_ht 20 ind_ht],...
   'Value',0,...
   'Callback','set_indepvar');
ind_ht = 20; tot_ht = tot_ht - ind_ht;
uicontrol('Parent',hiamrfig,...
   'Style','text',...
   'Position',[0 tot_ht 200 ind_ht],...
   'String', 'Regress on IAData using ILD_dist');
hregressia = uicontrol('Parent',hiamrfig,...
   'Style','checkbox',...
   'Position',[210 tot_ht 20 ind_ht],...
   'Value',1,...
   'Callback','set_indepvar');

%Commence analysis
ind_ht = 20; tot_ht = tot_ht - ind_ht - 30;
uicontrol('Parent',hiamrfig,...
   'Style','pushbutton',...
   'Position',[0 tot_ht 100 ind_ht],...
   'String', 'Commence Analysis',...
   'Callback','Engage_IAMultRegress');

%Cancel
uicontrol('Parent',hiamrfig,...
   'Style','pushbutton',...
   'Position',[iamrfig_wd-100 tot_ht 100 ind_ht],...
   'String', 'Cancel',...
   'Callback','close(hiamrfig)');

return