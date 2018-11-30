function varargout = SimpleSuicide_SD_fig_qi(varargin)
% SIMPLESUICIDE_SD_FIG_QI MATLAB code for SimpleSuicide_SD_fig_qi.fig
%      SIMPLESUICIDE_SD_FIG_QI, by itself, creates a new SIMPLESUICIDE_SD_FIG_QI or raises the existing
%      singleton*.
%
%      H = SIMPLESUICIDE_SD_FIG_QI returns the handle to a new SIMPLESUICIDE_SD_FIG_QI or the handle to
%      the existing singleton*.
%
%      SIMPLESUICIDE_SD_FIG_QI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMPLESUICIDE_SD_FIG_QI.M with the given input arguments.
%
%      SIMPLESUICIDE_SD_FIG_QI('Property','Value',...) creates a new SIMPLESUICIDE_SD_FIG_QI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SimpleSuicide_SD_fig_qi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SimpleSuicide_SD_fig_qi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SimpleSuicide_SD_fig_qi

% Last Modified by GUIDE v2.5 29-Nov-2018 17:06:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SimpleSuicide_SD_fig_qi_OpeningFcn, ...
                   'gui_OutputFcn',  @SimpleSuicide_SD_fig_qi_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SimpleSuicide_SD_fig_qi is made visible.
function SimpleSuicide_SD_fig_qi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SimpleSuicide_SD_fig_qi (see VARARGIN)

% Choose default command line output for SimpleSuicide_SD_fig_qi
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SimpleSuicide_SD_fig_qi wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SimpleSuicide_SD_fig_qi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function ar_Callback(hObject, eventdata, handles)
SimpleSuicide_SD_qi(handles)

function ac_Callback(hObject, eventdata, handles)
SimpleSuicide_SD_qi(handles)

function aa_Callback(hObject, eventdata, handles)
SimpleSuicide_SD_qi(handles)

function b0_Callback(hObject, eventdata, handles)
SimpleSuicide_SD_qi(handles)

function Ir_Callback(hObject, eventdata, handles)
SimpleSuicide_SD_qi(handles)

function Ic_Callback(hObject, eventdata, handles)
SimpleSuicide_SD_qi(handles)

function Ia_Callback(hObject, eventdata, handles)
SimpleSuicide_SD_qi(handles)

function Pn_Callback(hObject, eventdata, handles)
SimpleSuicide_SD_qi(handles)
