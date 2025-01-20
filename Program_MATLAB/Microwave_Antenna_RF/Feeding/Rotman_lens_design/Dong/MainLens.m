function varargout = MainLens(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton', gui_Singleton, ...
    'gui_OpeningFcn', @MainLens_OpeningFcn, ...
    'gui_OutputFcn', @MainLens_OutputFcn, ...
    'gui_LayoutFcn', [], ...
    'gui_Callback', []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State,varargin{:});
else
    gui_mainfcn(gui_State,varargin{:});
end

function MainLens_OpeningFcn(hObject,eventdata,handles,varargin)
handles.output = hObject;
handles.CenterF = 2;
set(handles.Center_Freq,'String',handles.CenterF);
handles.LoF = 1;        set(handles.Lower_Freq,'String',handles.LoF);
handles.UpF = 4;        set(handles.Upper_Freq,'String',handles.UpF);
handles.Eps = 1;        set(handles.Epso,'String',handles.Eps);
handles.thick = 0.001;  set(handles.Dthick,'String',handles.thick);
handles.Condu = 5.8e7;  set(handles.Conduc,'String',handles.Condu);
handles.Lstang = 0.0001;set(handles.Lstan,'String',handles.Lstang);
handles.TrmImp = 50;    set(handles.TermImp,'String',handles.TrmImp);
handles.Dspace = .5;    set(handles.Dspacing,'String',handles.Dspace);
handles.LenGama = 1;    set(handles.LensGama,'String',handles.LenGama);
handles.FocalAlph = 20; set(handles.FocalAlpha,'String',handles.FocalAlph);
handles.g3foc = 1.055;  set(handles.g3focal,'String',handles.g3foc);
handles.MaxThet = 30;   set(handles.MaxTheta,'String',handles.MaxThet);
handles.r3Len = 1;      set(handles.r3Lens,'String',handles.r3Len);
handles.eCentr = 0;     set(handles.eCentri,'String',handles.eCentr);
handles.MidLineLengt = 0;
set(handles.MidLineLength,'String',handles.MidLineLengt);
handles.TxEpss = 1;     set(handles.TxEps,'String',handles.TxEpss);
handles.Nbeams = 4;     set(handles.Nbeam,'String',handles.Nbeams);
handles.LBtapers = 1;   set(handles.LBtaper,'String',handles.LBtapers);
handles.Ndums = 1;      set(handles.Ndum,'String',handles.Ndums);
handles.BeamSiz = 5;
set(handles.BeamSize,'String',handles.BeamSiz);set(handles,BeamSizeSlider,'Value',handles.BeamSiz);
handles.DumSiz = 5;
set(handles.DumSize,'String',handles.DumSiz);set(handles.DumSizeSlider,'Value',handles.DumSiz);
handles.BeamNormals = 0;
set(handles.BeamNormal,'Value',0);set(handles.BeamOrigin,'Value',1);
handles.RcvNormals = 0;
set(handles.RcvNormal,'Value',0);set(handles.RcvMid,'Value',1);
handles.Nrcvs = 9; set(handles.Nrcv,'String',handles.Nrcvs);
handles.LRtapers = 1; set(handles.LRtaper,'String',handles.LRtapers);
handles.RrcvSiz = 5;
set(handles.RrcvSize,'String',handles.RrcvSiz);set(handles.RrcvSizeSlider,'Value',handles.RrcvSiz);
handles.Nsides = 10; set(handles.Nside,'String',handles.Nsides);
handles.tCurvatur = 0.96; set(handles.tCurvature,'String',handles.tCurvatur);
handles.LDtapers = 1; set(handles.LDtaper,'String',handles.LDtapers);
handles.Focal4Alph = 1; set(handles.Focal4Alpha,'String',handles.Focal4Alph);
handles.Focal4Bet = 1; set(handles.Focal4Beta,'String',handles.Focal4Bet);
handles.Max4Thet = 40; set(handles.Max4Theta,'String',handles.Max4Thet);
handles.r4Len = 1; set(handles.r4Lens,'String',handles.r4Len);
handles.TransmDif = 'Tx Line length';
set(handles.TransmDiff,'String',handles.TransmDif);
handles.Lamda = 0.3/handles.CenterF; handles.WidthTrans = micstripW(handles.Eps,handles.thick,handles.TrmImp);
handles.Pinfo = 0;
handles.Xite = 0;
handles.Wtrans = 0;
setappdata(0,'hMainLens',gcf);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Dspace',handles.Dspace);
setappdata(hMainLens,'CenterF',handles.CenterF);
setappdata(hMainLens,'Wtrans',handles.Wtrans);
setappdata(hMainLens,'WidthTrans',handles.WidthTrans);
setappdata(hMainLens,'tCurvatur',handles.tCurvatur);
setappdata(hMainLens,'LBtapers',handles.LBtapers);
setappdata(hMainLens,'LDtapers',handles.LDtapers);
setappdata(hMainLens,'LRtapers',handles.LRtapers);
setappdata(hMainLens,'BeamSiz',handles.BeamSiz);
setappdata(hMainLens,'DumSiz',handles.DumSiz);
setappdata(hMainLens,'RrcvSiz',handles.RrcvSiz);
setappdata(hMainLens,'eCentr',handles.eCentr);
setappdata(hMainLens,'LenGama',handles.LenGama);
setappdata(hMainLens,'FocalAlph',handles.FocalAlph);
setappdata(hMainLens,'g3foc',handles.g3foc);
setappdata(hMainLens,'MaxThet',handles.MaxThet);
setappdata(hMainLens,'r3Len',handles.r3Len);
setappdata(hMainLens,'Nbeams',handles.Nbeams);
setappdata(hMainLens,'Ndums',handles.Ndums);
setappdata(hMainLens,'Nsides',handles.Nsides);
setappdata(hMainLens,'Nrcvs',handles.Nrcvs);
setappdata(hMainLens,'BeamNormals',handles.BeamNormals);
setappdata(hMainLens,'RcvNormals',handles.RcvNormals);

setappdata(hMainLens,'Status','Welcome to Rotman Lens Designer!');
set(handles.Status,'String',getappdata(hMainLens,'Status'),'FontSize',18.0);

[handles.Pinfo,handles.Xite,handles.Wtrans] = strucLens0(handles.Dspace,handles.CenterF,handles.WidthTrans,handles.tCurvatur,handles.LBtapers,handles.LDtapers,handles.LRtapers,...
    handles.BeamSiz,handles.DumSiz,handles.RrcvSiz,handles.eCentr,handles.LenGama,handles.FocalAlph,handles.g3foc,handles.MaxThet,handles.r3Len,handles.Nbeams,handles.Ndums,...
    handles.Nsides,handles.Nrcvs,handles.BeamNormals,handles.RcvNormals);

setappdata(hMainLens,'WidthTrans',handles.WidthTrans);
setappdata(hMainLens,'Lamda',handles.Lamda);
setappdata(hMainLens,'Pinfo',handles.Pinfo);
setappdata(hMainLens,'Xite',handles.Xite);
setappdata(hMainLens,'Wtrans',handles.Wtrans);
setappdata(hMainLens,'NN',1);
setappdata(hMainLens,'SCcontrl',1);
setappdata(hMainLens,'TrmImp',handles.TrmImp);
setappdata(hMainLens,'Eps',handles.Eps);
setappdata(hMainLens,'thick',handles.thick);
setappdata(hMainLens,'Lstang',handles.Lstang);
setappdata(hMainLens,'SCFreq',handles.CenterF);
setappdata(hMainLens,'fhSCoupling',@SCoupling);
setappdata(hMainLens,'fhS2Coupling',@S2Coupling);
setappdata(hMainLens,'NN2',1);
setappdata(hMainLens,'MM2',1);
setappdata(hMainLens,'LowFreq2',handles.LoF);
setappdata(hMainLens,'HighFreq2',handles.UpF);
setappdata(hMainLens,'Nfre2',100);
setappdata(hMainLens,'CtrlRefe',1);
setappdata(hMainLens,'CtrlRefe2',1);
guidata(hObject,handles);

function varargout = MainLens_OutputFcn(hObject,eventdata,handles)
varargout{1} = handles.output;

function r4Lens_Callback(hObject,eventdata,handles)
handles.r4Len = str2double(get(handles.r4Lens,'String'));
set(handles.r4Lens,'String',handles.r4Len);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'r4len',handles.r4Len);
guidata(hObject,handles);

function r4Lens_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Upper_Freq_Callback(hObject,eventdata,handles)
handles.UpF = str2double(get(handles.Upper_Freq,'String'));
set(handles.Upper_Freq,'String',handles.UpF);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'UpF',handles.UpF);
guidata(hObject,handles);

function Upper_Freq_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nbeam_Callback(hObject,eventdata,handles)
handles.Nbeams = str2double(get(handles.Nbeam,'String'));
set(handles.Nbeam,'String',handles.Nbeams);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Nbeams',handles.Nbeams);
guidata(hObject,handles);

function Nbeam_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LensGama_Callback(hObject,eventdata,handles)
handles.LenGama = str2double(get(handles.LensGama,'String'));
set(handles.LensGama,'String',handles.LenGama);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'LenGama',handles.LenGama);
guidata(hObject,handles);

function LensGama_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function g3focal_Callback(hObject,eventdata,handles)
handles.g3foc = str2double(get(handles.g3focal,'String'));
set(handles.g3focal,'String',handles.g3foc);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'g3foc',handles.g3foc);
guidata(hObject,handles);

function g3focal_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaxTheta_Callback(hObject,eventdata,handles)
handles.MaxThet = str2double(get(handles.MaxTheta,'String'));
set(handles.MaxTheta,'String',handles.MaxThet);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'MaxThet',handles.MaxThet);
guidata(hObject,handles);

function MaxTheta_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Conduc_Callback(hObject,eventdata,handles)
handles.Condu = str2double(get(handles.Conduc,'String'));
set(handles.Conduc,'String',handles.Condu);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Condu',handles.Condu);
guidata(hObject,handles);

function Conduc_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Epso_Callback(hObject,eventdata,handles)
handles.Eps = str2double(get(handles.Epso,'String'));
set(handles.Epso,'String',handles.Eps);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Eps',handles.Eps);
guidata(hObject,handles);

function Epso_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Dthick_Callback(hObject,eventdata,handles)
handles.thick = str2double(get(handles.Dthick,'String'));
set(handles.Dthick,'String',handles.thick);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'thick',handles.thick);
guidata(hObject,handles);

function Dthick_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nside_Callback(hObject,eventdata,handles)
handles.Nsides = str2double(get(handles.Nside,'String'));
set(handles.Nside,'String',handles.Nsides);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Nsides',handles.Nsides);
guidata(hObject,handles);

function Nside_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tCurvature_Callback(hObject,eventdata,handles)
handles.tCurvatur = str2double(get(handles.tCurvature,'String'));
set(handles.tCurvature,'String',handles.tCurvatur);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'tCurvatur',handles.tCurvatur);
guidata(hObject,handles);

function tCurvature_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r3Lens_Callback(hObject,eventdata,handles)
handles.r3Len = str2double(get(handles.r3Lens,'String'));
set(handles.r3Lens,'String',handles.r3Len);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'r3Len',handles.r3Len);
guidata(hObject,handles);

function r3Lens_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nrcv_Callback(hObject,eventdata,handles)
handles.Nrcvs = str2double(get(handles.Nrcv,'String'));
set(handles.Nrcv,'String',handles.Nrcvs);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Nrcvs',handles.Nrcvs);
guidata(hObject,handles);

function Nrcv_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Center_Freq_Callback(hObject,eventdata,handles)
handles.CenterF = str2double(get(handles.Center_Freq,'String'));
set(handles.Center_Freq,'String',handles.CenterF);

hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'CenterF',handles.CenterF);
handles.WidthTrans = micstripW(handles.Eps,handles.thick,handles.TrmImp);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'WidthTrans',handles.WidthTrans);
guidata(hObject,handles);

function Lower_Freq_Callback(hObject,eventdata,handles)
handles.LoF = str2double(get(handles.Lower_Freq,'String'));
set(handles.Lower_Freq,'String',handles.LoF);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'LoF',handles.LoF);
guidata(hObject,handles);

function FSave_Callback(hObject,eventdata,handles)
[filename,pathname] = uiputfile;
save(filename,'-struct','handles');
guidata(hObject,handles);

function ToolSpara_Callback(hObject,eventdata,handles)

MenuSpara

function [Snew] = SCoupling
hMainLens = getappdata(0,'hMainLens');
Wtrans    = getappdata(hMainLens,'Wtrans');
SCFreq    = getappdata(hMainLens,'SCFreq');
CtrlRefe  = getappdata(hMainLens,'CtralRefe');

NN        = getappdata(hMainLens,'NN');
TrmImp    = getappdata(hMainLens,'TrmImp');
Eps       = getappdata(hMainLens,'Eps');
thick     = getappdata(hMainLens,'thick');
Lstang    = getappdata(hMainLens,'Lstang');
SCcontrl = getappdata(hMainLens,'SCcontrl');
Pinfo     = getappdata(hMainLens,'Pinfo');
Xite      = getappdata(hMainLens,'Xite');

[Snew] = Menu_Coupling(NN,Wtrans,TrmImp,Eps,thick,Lstang,SCFreq,Pinfo,Xite,SCcontrl,CtrlRefe);

function S2Coupling
hMainLens    = getappdata(0,'hMainLens');
Wtrans       = getappdata(hMainLens,'Wtrans');
TrmImp       = getappdata(hMainLens,'TrmImp');
Eps          = getappdata(hMainLens,'Eps');
thick        = getappdata(hMainLens,'thick');
Lstang       = getappdata(hMainLens,'Lstang');
SCcontrl2   = getappdata(hMainLens,'SCcontrl2');
Pinfo        = getappdata(hMainLens,'Pinfo');
Xite         = getappdata(hMainLens,'Xite');
NN2          = getappdata(hMainLens,'NN2');
MM2          = getappdata(hMainLens,'MM2');
LowFreq2     = getappdata(hMainLens,'LowFreq2');
HighFreq2    = getappdata(hMainLens,'HighFreq2');
Nfre2        = getappdata(hMainLens,'Nfre2');
CtrlRefe2    = getappdata(hMainLens,'CtrlRefe2');
Menu_Coupling2(NN2,MM2,Wtrans,LowFreq2,HighFreq2,Nfre2,TrmImp,Eps,thick,Lstang,Pinfo,Xite,SCcontrl2,CtrlRefe2);

function FocalAlpha_Callback(hObject,eventdata,handles)
handles.FocalAlph = str2double(get(handles.FocalAlpha,'String'));
set(handles.FocalAlpha,'String',handles.FocalAlph);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'FocalAlph',handles.FocalAlph);
guidata(hObject,handles);

function FocalAlpha_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Focal4Beta_Callback(hObject,eventdata,handles)
handles.Focal4Bet = str2double(get(handles.Focal4Beta,'String'));
set(handles.Focal4Beta,'String',handles.Focal4Bet);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Focal4Bet',handles.Focal4Bet);
guidata(hObject,handles);

function Focal4Beta_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Max4Theta_Callback(hObject,eventdata,handles)
handles.Max4Thet = str2double(get(handles.Max4Theta,'String'));
set(handles.Max4Theta,'String',handles.Max4Thet);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Max4Thet',handles.Max4Thet);
guidata(hObject,handles);

function Max4Theta_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function eCentri_Callback(hObject,eventdata,handles)
handles.eCentr = str2double(get(handles.eCentri,'String'));
set(handles.eCentri,'String',handles.eCentr);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'eCentr',handles.eCentr);
guidata(hObject,handles);

function eCentri_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BeamNormal_Callback(hObject,eventdata,handles)
if (get(handles.BeamNormal,'Value') == get(handles.BeamNormal,'Max'))
    set(handles.BeamOrigin,'Value',0);
    handles.BeamNormals = 1;
else
    set(handles.BeamOrigin,'Value',1);
    handles.BeamNormals = 0;
end
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'BeamNormals',handles.BeamNormals);
guidata(hObject,handles);

function BeamOrigin_Callback(hObject,eventdata,handles)
if (get(handles.BeamOrigin,'Value') == get(handles.BeamOrigin,'Max'))
    set(handles.BeamNormal,'Value',0);
    handles.BeamNormals = 0;
else
    set(handles.BeamNormal,'Value',1);
    handles.BeamNormals = 1;
end
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'BeamNormals',handles.BeamNormals);
guidata(hObject,handles);

function RcvNormal_Callback(hObject,eventdata,handles)
if (get(handles.RcvNormal,'Value') == get(handles.RcvNormal,'Max'))
    set(handles.RcvMid,'Value',0);
    handles.RcvNormals = 0;
else
    set(handles.RcvMid,'Value',1);
    handles.RcvNormals = 1;
end
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'RcvNormals',handles.RcvNormals);
guidata(hObject,handles);

function RcvMid_Callback(hObject,eventdata,handles)
if (get(handles.RcvMid,'Value') == get(handles.RcvMid,'Max'))
    set(handles.RcvNormal,'Value',0);
    handles.RcvNormals = 1;
else
    set(handles.RcvNormal,'Value',1);
    handles.RcvNormals = 1;
end
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'RcvNormals',handles.RcvNormals);
guidata(hObject,handles);

function LRtaper_Callback(hObject,eventdata,handles)
handles.LRtapers = str2double(get(handles.LRtaper,'String'));
set(handles.LRtaper,'String',handles.LRtapers);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'LRtapers',handles.LRtapers);
guidata(hObject,handles);

function LRtaper_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ndum_Callback(hObject,eventdata,handles)
handles.Ndums = str2double(get(handles.Ndum,'String'));
set(handles.Ndum,'String',handles.Ndums);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Ndums',handles.Ndums);
guidata(hObject,handles);

function Ndum_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BeamSize_Callback(hObject,eventdata,handles)
handles.BeamSiz = round(str2double(get(handles.BeamSize,'String')));
handles.DumSiz  = round(str2double(get(handles.DumSize,'String')));
set(handles.BeamSize,'String',handles.BeamSiz);
set(handles.BeamSizeSlider,'Value',handles.BeamSiz);
if handles.BeamSiz+handles.DumSiz>10
    handles.DumSiz = 10-handles.BeamSiz;
    set(handles.DumSize,'String',handles.DumSiz);
    set(handles.DumSizeSlider,'Value',handles.DumSiz);
end
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'BeamSiz',handles.BeamSiz);
setappdata(hMainLens,'DumSiz',handles.DumSiz);
guidata(hObject,handles);

function BeamSize_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BeamSizeSlider_Callback(hObject,eventdata,handles)
handles.BeamSiz = round(get(handles.BeamSizeSlider,'Value'));
handles.DumSiz  = round(str2double(get(handles.DumSize,'String')));
set(handles.BeamSize,'String',handles.BeamSiz);

if handles.BeamSiz+handles.DumSiz>10
    handles.DumSiz = 10-handles.BeamSiz;
    set(handles.DumSize,'String',handles.DumSiz);
    set(handles.DumSizeSlider,'Value',handles.DumSiz);
end
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'BeamSiz',handles.BeamSiz);
setappdata(hMainLens,'DumSiz',handles.DumSiz);

function BeamSizeSlider_CreateFcn(hObject,eventdata,handles)
if  isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function RrcvSize_Callback(hObject,eventdata,handles)
handles.RrcvSiz = round(str2double(get(handles.RrcvSize,'String')));
set(handles.RrcvSize,'String',handles.RrcvSiz);
set(handles.RrcvSizeSlider,'Value',handles.RrcvSiz);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'RrcvSiz',handles.RrcvSiz);
guidata(hObject,handles);

function RrcvSizeSlider_Callback(hObject,eventdata,handles)
handles.RrcvSiz = round(get(handles.RrcvSizeSlider,'Value'));
set(handles.RrcvSize,'String',handles.RrcvSiz);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'RrcvSiz',handles.RrcvSiz);
guidata(hObject,handles);

function RrcvSize_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Lstan_Callback(hObject,eventdata,handles)
handles.Lstang = str2double(get(handles.Lstan,'String'));
set(handles.Lstan,'String',handles.Lstang);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Lstang',handles.Lstang);
guidata(hObject,handles);

function Lstan_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TermImp_Callback(hObject,eventdata,handles)
handles.TrmImp = str2double(get(handles.TermImp,'String'));
set(handles.TermImp,'String',handles.TrmImp);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'TrmImp',handles.TrmImp);
guidata(hObject,handles);

function Dspacing_Callback(hObject,eventdata,handles)
handles.Dspace = str2double(get(handles.Dspacing,'String'));
set(handles.Dspacing,'String',handles.Dspace);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Dspace',handles.Dspace);
guidata(hObject,handles);

function LBtaper_Callback(hObject,eventdata,handles)
handles.LBtapers = str2double(get(handles.LBtaper,'String'));
set(handles.LBtaper,'String',handles.LBtapers);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'LBtapers',handles.LBtapers);
guidata(hObject,handles);

function Focal4Alpha_Callback(hObject,eventdata,handles)
handles.Focal4Alph = str2double(get(handles.Focal4Alpha,'String'));
set(handles.Focal4Alpha,'String',handles.Focal4Alph);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'Focal4Alph',handles.Alph);
guidata(hObject,handles);

function DumSize_Callback(hObject,eventdata,handles)
handles.DumSiz = round(str2double(get(handles.DumSize,'String')));
handles.BeamSiz = round(str2double(get(handles.BeamSize,'String')));
set(handles.DumSize,'String',handles.DumSiz);
set(handles.DumSizeSlider,'Value',handles.DumSiz);
if handles.BeamSiz+handles.DumSiz>10
    handles.BeamSiz = 10-handles.DumSiz;
    set(handles.BeamSize,'String',handles.BeamSiz);
    set(handles.BeamSizeSlider,'Value',handles.BeamSiz);
end
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'DumSiz',handles.DumSiz);
setappdata(hMainLens,'BeamSiz',handles.BeamSiz);
guidata(hObject,handles);

function DumSize_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DumSizeSlider_Callback(hObject,eventdata,handles)
handles.DumSiz = round(get(handles.DumSizeSlider,'Value'));
handles.BeamSiz = round(str2double(get(handles.BeamSize,'String')));
set(handles.DumSize,'String',handles.DumSiz);
if handles.BeamSiz+handles.DumSiz>10
    handles.BeamSiz = 10-handles.DumSiz;
    set(handles.BeamSize,'String',handles.BeamSiz);
    set(handles.BeamSizeSlider,'Value',handles.BeamSiz);
end
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'DumSiz',handles.DumSiz);
setappdata(hMainLens,'BeamSiz',handles.BeamSiz);
guidata(hObject,handles);

function DumSizeSlider_CreateFcn(hObject,eventdata,handles)
if isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function RefreshButton_Callback(hObject,eventdata,handles)
handles.WidthTrans = micstripW(handles.Eps,handles.thick,handles.TrmImp);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'WidthTrans',handles.WidthTrans);
off;
hMainLens = getappdata(0,'hMainLens');
Dspace = getappdata(hMainLens,'Dspace');
CenterF = getappdata(hMainLens,'CenterF');
WidthTrans = getappdata(hMainLens,'WidthTrans');
tCurvatur = getappdata(hMainLens,'tCurvatur');
LBtapers = getappdata(hMainLens,'LBtapers');
LDtapers = getappdata(hMainLens,'LDtapers');
LRtapers = getappdata(hMainLens,'LRtapers');
BeamSiz = getappdata(hMainLens,'BeamSiz');
DumSiz = getappdata(hMainLens,'DumSiz');
RrcvSiz = getappdata(hMainLens,'RrcvSiz');
eCentr = getappdata(hMainLens,'eCentr');
LenGama = getappdata(hMainLens,'LenGama');
FocalAlph = getappdata(hMainLens,'FocalAlph');
g3foc = getappdata(hMainLens,'g3foc');
MaxThet = getappdata(hMainLens,'MaxThet');
r3Len = getappdata(hMainLens,'r3Len');
Nbeams = getappdata(hMainLens,'Nbeams');
Ndums = getappdata(hMainLens,'Ndums');
Nsides = getappdata(hMainLens,'Nsides');
Nrcvs = getappdata(hMainLens,'Nrcvs');
BeamNormals = getappdata(hMainLens,'BeamNormals');
RcvNormals = getappdata(hMainLens,'RcvNormals');

[handles.Pinfo,handles.Xite,handles.Wtrans] = strucLens(Dspace,CenterF,WidthTrans,tCurvatur,LBtapers,LDtapers,LRtapers,BeamSiz,DumSiz,RrcvSiz,eCentr,LenGama,FocalAlph,...
    g3foc,MaxThet,r3Len,Nbeams,Ndums,Nsides,Nrcvs,BeamNormals,RcvNormals);
setappdata(hMainLens,'Pinfo',handles.Pinfo);
setappdata(hMainLens,'Xite',handles.Xite);

axis equal;
handles.TransmDif = handles.MidLineLengt+handles.Wtrans;
set(handles.TransmDiff,'String',handles.TransmDif);
setappdata(hMainLens,'Wtrans',handles.TransmDif);
setappdata(hMainLens,'Status','Please click Refresh button after you update the lens parameters.');
set(handles.Status,'String',getappdata(hMainLens,'Status'),'FontSize',15.0);
guidata(hObject,handles);

function LDtaper_Callback(hObject,eventdata,handles)
handles.LDtapers = str2double(get(handles.LDtaper,'String'));
set(handles.LDtaper,'String',handles.LDtapers);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'LDtapers',handles.LDtapers);
guidata(hObject,handles);

function LDtaper_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TransmDiff_Callback(hObject,eventdata,handles)
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'TransmDiffs',handles.TransmDiffs);

function TransmDiff_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MidLineLength_Callback(hObject,eventdata,handles)
handles.MidLineLengt = str2double(get(handles.MidLineLength,'String'));
set(handles.MidLineLength,'String',handles.MidLineLengt);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'MidLineLengt',handles.MidLineLengt);
guidata(hObject,handles);

function MidLineLength_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TxEps_Callback(hObject,eventdata,handles)
handles.TxEpss = str2double(get(handles.TxEps,'String'));
set(handles.TxEps,'String',handles.TxEpss);
hMainLens = getappdata(0,'hMainLens');
setappdata(hMainLens,'TxEpss',handles.TxEpss);
guidata(hObject,handles);

function TxEps_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end