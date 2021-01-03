function varargout = EquationSolver(varargin)
% EQUATIONSOLVER MATLAB code for EquationSolver.fig
%      EQUATIONSOLVER, by itself, creates a new EQUATIONSOLVER or raises the existing
%      singleton*.
%
%      H = EQUATIONSOLVER returns the handle to a new EQUATIONSOLVER or the handle to
%      the existing singleton*.
%
%      EQUATIONSOLVER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EQUATIONSOLVER.M with the given input arguments.
%
%      EQUATIONSOLVER('Property','Value',...) creates a new EQUATIONSOLVER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EquationSolver_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EquationSolver_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EquationSolver

% Last Modified by GUIDE v2.5 03-Jan-2021 16:31:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EquationSolver_OpeningFcn, ...
                   'gui_OutputFcn',  @EquationSolver_OutputFcn, ...
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


% --- Executes just before EquationSolver is made visible.
function EquationSolver_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EquationSolver (see VARARGIN)

% Choose default command line output for EquationSolver
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EquationSolver wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EquationSolver_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function fileName_Callback(hObject, eventdata, handles)
% hObject    handle to fileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileName as text
%        str2double(get(hObject,'String')) returns contents of fileName as a double


% --- Executes during object creation, after setting all properties.
function fileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in readFile.
function readFile_Callback(hObject, eventdata, handles)
set(handles.readFile, 'value', 1);
set(handles.error,'string',"");
if get(handles.writeEquations, 'value') == 1
   set(handles.writeEquations, 'value', 0); 
   set(handles.equations, 'enable', 'off');
   set(handles.fileName, 'enable' , 'on');
   set(handles.methods, 'enable', 'off');
end


% --- Executes on button press in solve.
function solve_Callback(hObject, eventdata, handles)
parser = EquationParser;
solver = LinearSolver;
if get(handles.readFile, 'value') == 1
    if get(handles.fileName, 'string') == "" || exist(get(handles.fileName, 'string'), 'file') ~= 2
        set(handles.error, 'string', "enter a correct file name");
        return;
    end
    fid = fopen(get(handles.fileName, 'string'));
    txt = textscan(fid, "%s", 'delimiter', '\n');
    num_eqns = txt{1}{1};
    method_name = txt{1}{2};
    if strcmpi(method_name, "Gaussian-elimination")
        method_value = 1;
    elseif strcmpi(method_name, "Gaussian-Jordan")
        method_value = 2;
    elseif strcmpi(method_name, "LU decomposition")
        method_value = 3;
    elseif strcmpi(method_name, "Gauss-Seidel")
        method_value = 4;
    else
        set(handles.error,'string',"wrong file format");
        return;
    end
    set(handles.error,'string',"");
    j = 4;
    s = char(txt{1}{3});
    for i=1:str2double(num_eqns)-1
        s = char(s, txt{1}{j});
        j = j + 1;
    end
    if method_value == 4
        str22 = txt{1}{j};
        str22 = deblank(str22);
        intial_guess = parser.getInitialGuesses(str22);
    end
    fclose(fid);
else
s = get(handles.equations, 'string');
if isempty(deblank(s))
    set(handles.error, 'string' ,'please enter the equations');
    return;
end
set(handles.error, 'string' ,'');
method_value = get(handles.methods, 'value');
if method_value == 4
    str22 = get(handles.initial_guesses, 'string');
    if str22 == ""
        set(handles.error, 'string' ,'please enter intial guesses');
        return;
    end
    set(handles.error, 'string' ,'');
    intial_guess = parser.getInitialGuesses(str22);
end
% [A,b, index_to_var] = parser.equationsToMatrix(s);
end
mi = get(handles.max_iter, 'string');
if mi == ""
    max_iter = 50;
else
    max_iter = str2double(mi);
end
ee = get(handles.eps, 'string');
if ee == ""
    eps = 0.00001;
else
    eps = str2double(ee);
end
[A,b, index_to_var] = parser.equationsToMatrix(s);
tic
% [A,b, index_to_var] = parser.equationsToMatrix(s);
switch method_value
    case 1
        [x, root_str] = solver.Gauss(A, b, eps, max_iter, index_to_var);
        set(handles.numIterations, 'string' , "");
        set(handles.ea, 'string' , "");
        set(handles.all_iterations, 'string', "");
    case 2
        [x, root_str] = solver.Gauss_Jordan(A, b, index_to_var);
        set(handles.numIterations, 'string' , "");
        set(handles.ea, 'string' , "");
        set(handles.all_iterations, 'string', "");
    case 3
        [x, root_str] = solver.LU_Decomp(A, b, index_to_var);
        set(handles.numIterations, 'string' , "");
        set(handles.ea, 'string' , "");
        set(handles.all_iterations, 'string', "");
    case 4
        [x, iter_str, num_iter, err, root_str] = solver.Gauss_Seidel(A, b, eps, max_iter, intial_guess, index_to_var);
        set(handles.all_iterations, 'string', iter_str);
        if num_iter >= max_iter
            set(handles.numIterations, 'string' , "DID NOT CONVERGE!!!");
        else
            set(handles.numIterations, 'string' , sprintf("number of iterations = %d",num_iter));    
        end
        set(handles.numIterations, 'string' , sprintf("number of iterations = %d",num_iter));
        set(handles.ea, 'string' , sprintf("finished with error = %f", err));
end
for i = 1:length(x)
    root_str = strcat(root_str, sprintf("%f         ", x(i))); 
end
set(handles.solvedRoots, 'string' , root_str);
time = toc;
set(handles.timeTaken, 'string' , sprintf("Time taken = %f", time));
fileID = fopen("output.txt", 'w');
fprintf(fileID,"roots:\n%s\n", root_str);
if method_value == 4
   fprintf(fileID,"All iterations:\n%s",iter_str);
   if num_iter >= max_iter
        fprintf(fileID,"DID NOT CONVERGE!!!!\n");
   else
        fprintf(fileID,"number of Iterations = %d\n",num_iter);       
   end

   fprintf(fileID, "finished with error = %f\n",err);
end
fprintf(fileID,"time taken = %f\n", time);
fclose(fileID);






% --- Executes on button press in writeEquations.
function writeEquations_Callback(hObject, eventdata, handles)
set(handles.writeEquations, 'value', 1);
if get(handles.readFile, 'value') == 1
    set(handles.readFile, 'value', 0);
    set(handles.error,'string',"");
    set(handles.equations, 'enable', 'on');
    set(handles.fileName, 'enable', 'off');
    set(handles.methods, 'enable', 'on');
end



function equations_Callback(hObject, eventdata, handles)
% hObject    handle to equations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of equations as text
%        str2double(get(hObject,'String')) returns contents of equations as a double


% --- Executes during object creation, after setting all properties.
function equations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to equations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in methods.
function methods_Callback(hObject, eventdata, handles)
if get(handles.methods, 'value') == 4
   set(handles.initial_guesses, 'enable', 'on');
else
    set(handles.initial_guesses, 'enable', 'off'); 
end


% --- Executes during object creation, after setting all properties.
function methods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initial_guesses_Callback(hObject, eventdata, handles)
% hObject    handle to initial_guesses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initial_guesses as text
%        str2double(get(hObject,'String')) returns contents of initial_guesses as a double


% --- Executes during object creation, after setting all properties.
function initial_guesses_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initial_guesses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eps_Callback(hObject, eventdata, handles)
% hObject    handle to eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eps as text
%        str2double(get(hObject,'String')) returns contents of eps as a double


% --- Executes during object creation, after setting all properties.
function eps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_iter_Callback(hObject, eventdata, handles)
% hObject    handle to max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_iter as text
%        str2double(get(hObject,'String')) returns contents of max_iter as a double


% --- Executes during object creation, after setting all properties.
function max_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in all_iterations.
function all_iterations_Callback(hObject, eventdata, handles)
% hObject    handle to all_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns all_iterations contents as cell array
%        contents{get(hObject,'Value')} returns selected item from all_iterations


% --- Executes during object creation, after setting all properties.
function all_iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to all_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in solvedRoots.
function solvedRoots_Callback(hObject, eventdata, handles)
% hObject    handle to solvedRoots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns solvedRoots contents as cell array
%        contents{get(hObject,'Value')} returns selected item from solvedRoots


% --- Executes during object creation, after setting all properties.
function solvedRoots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to solvedRoots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
