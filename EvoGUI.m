function varargout = EvoGUI(varargin)
    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @EvoGUI_OpeningFcn, ...
                       'gui_OutputFcn',  @EvoGUI_OutputFcn, ...
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


% --- Executes just before EvoGUI is made visible.
function EvoGUI_OpeningFcn(hObject, eventdata, handles, varargin)
    % Choose default command line output for EvoGUI
    handles.output = hObject;

    sampleFiles = cell(0,0); %Pass a blank terrain sample array
    handles.evoTerr = EvoClass(sampleFiles,513,5,0.5,'spline',8,0.5,0.1);
    handles.evoTerr.initEvo();

    handles.genCounter = 0;
    handles.axesArray = [handles.axes1 handles.axes2 handles.axes3 handles.axes4 handles.axes5 handles.axes6 handles.axes7 handles.axes8];
    handles.checkboxArray = [handles.checkbox1 handles.checkbox2 handles.checkbox3 handles.checkbox4 handles.checkbox5 handles.checkbox6 handles.checkbox7 handles.checkbox8];
    for iter = 1:8
        surf(handles.axesArray(iter), handles.evoTerr.getHeightMap(iter), 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
        caxis(handles.axesArray(iter),[-0.01 1]);
        cm = [0.839 0.714 0.605]; %pink
        colormap(handles.axesArray(iter), cm);
        %view(handles.axesArray(iter), 140, 30); % Matlab renders upside down, rotate the view to be inline with matrix orientation 
        light('Parent', handles.axesArray(iter)); %specifically create a light for each object
        %camlight(lh,-45,30);
        %if iter == 1 %Stops a bug that causes the last graph to be brighter then the rest
            %camlight(lh, 0,20); 
        %end
        camproj(handles.axesArray(iter), 'perspective');
        camzoom(handles.axesArray(iter), 1.8);
        daspect(handles.axesArray(iter), [125 125 1])
        axis(handles.axesArray(iter), 'off');
    end

    %rotate3d on;

    % Update handles structure
    guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = EvoGUI_OutputFcn(hObject, eventdata, handles) 
    % Get default command line output from handles structure
    varargout{1} = handles.output;


function checkbox1_Callback(hObject, eventdata, handles)
    if get(hObject, 'Value') == get(hObject, 'Max')
        %Check box is ticked 
        %setAsParent will return false if the parentLimit is reached
        if ~(handles.evoTerr.setAsParent(1,true))
           set(hObject, 'Value', get(hObject, 'Min')); 
        end
    else
        %Check box is unticked
        handles.evoTerr.setAsParent(1,false);        
    end

    % Update handles structure
    guidata(hObject, handles);

function checkbox2_Callback(hObject, eventdata, handles)
    if get(hObject, 'Value') == get(hObject, 'Max')
        %Check box is ticked 
        %setAsParent will return false if the parentLimit is reached
        if ~(handles.evoTerr.setAsParent(2,true))
           set(hObject, 'Value', get(hObject, 'Min')); 
        end
    else
        %Check box is unticked
        handles.evoTerr.setAsParent(2,false);        
    end

    % Update handles structure
    guidata(hObject, handles);


function checkbox3_Callback(hObject, eventdata, handles)
    if get(hObject, 'Value') == get(hObject, 'Max')
        %Check box is ticked 
        %setAsParent will return false if the parentLimit is reached
        if ~(handles.evoTerr.setAsParent(3,true))
           set(hObject, 'Value', get(hObject, 'Min')); 
        end
    else
        %Check box is unticked
        handles.evoTerr.setAsParent(3,false);        
    end

    % Update handles structure
    guidata(hObject, handles);


function checkbox4_Callback(hObject, eventdata, handles)
    if get(hObject, 'Value') == get(hObject, 'Max')
        %Check box is ticked 
        %setAsParent will return false if the parentLimit is reached
        if ~(handles.evoTerr.setAsParent(4,true))
           set(hObject, 'Value', get(hObject, 'Min')); 
        end
    else
        %Check box is unticked
        handles.evoTerr.setAsParent(4,false);        
    end

    % Update handles structure
    guidata(hObject, handles);


function checkbox5_Callback(hObject, eventdata, handles)
    if get(hObject, 'Value') == get(hObject, 'Max')
        %Check box is ticked 
        %setAsParent will return false if the parentLimit is reached
        if ~(handles.evoTerr.setAsParent(5,true))
           set(hObject, 'Value', get(hObject, 'Min')); 
        end
    else
        %Check box is unticked
        handles.evoTerr.setAsParent(5,false);        
    end

    % Update handles structure
    guidata(hObject, handles);


function checkbox6_Callback(hObject, eventdata, handles)
    if get(hObject, 'Value') == get(hObject, 'Max')
        %Check box is ticked 
        %setAsParent will return false if the parentLimit is reached
        if ~(handles.evoTerr.setAsParent(6,true))
           set(hObject, 'Value', get(hObject, 'Min')); 
        end
    else
        %Check box is unticked
        handles.evoTerr.setAsParent(6,false);        
    end

    % Update handles structure
    guidata(hObject, handles);

function checkbox7_Callback(hObject, eventdata, handles)
    if get(hObject, 'Value') == get(hObject, 'Max')
        %Check box is ticked 
        %setAsParent will return false if the parentLimit is reached
        if ~(handles.evoTerr.setAsParent(7,true))
           set(hObject, 'Value', get(hObject, 'Min')); 
        end
    else
        %Check box is unticked
        handles.evoTerr.setAsParent(7,false);        
    end

    % Update handles structure
    guidata(hObject, handles);

function checkbox8_Callback(hObject, eventdata, handles)
    if get(hObject, 'Value') == get(hObject, 'Max')
        %Check box is ticked 
        %setAsParent will return false if the parentLimit is reached
        if ~(handles.evoTerr.setAsParent(8,true))
           set(hObject, 'Value', get(hObject, 'Min')); 
        end
    else
        %Check box is unticked
        handles.evoTerr.setAsParent(8,false);        
    end

    % Update handles structure
    guidata(hObject, handles);

function evoButton_Callback(hObject, eventdata, handles)
    handles.evoTerr.evolve();

    handles.genCounter = handles.genCounter + 1;
    strGen = num2str(handles.genCounter);
    set(handles.textGeneration,'String',strGen);
    
    for iter = 1:8
        surf(handles.axesArray(iter), handles.evoTerr.getHeightMap(iter), 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
        caxis(handles.axesArray(iter),[0 1]);
        cm = [0.839 0.714 0.605]; %pink
        colormap(handles.axesArray(iter), cm);
        %view(handles.axesArray(iter), 140, 30);
        lh = light('Parent', handles.axesArray(iter)); %specifically create a light for each object
        %camlight(lh, 0,20);
        camproj(handles.axesArray(iter), 'perspective');
        camzoom(handles.axesArray(iter), 1.8);
        daspect(handles.axesArray(iter), [125 125 1])
        axis(handles.axesArray(iter), 'off');
    end
    
    % Reset the checkboxes
    for iter = 1:8
        if get(handles.checkboxArray(iter), 'Value') == get(handles.checkboxArray(iter), 'Max');
            set(handles.checkboxArray(iter), 'Value', get(handles.checkboxArray(iter), 'Min'));
            handles.evoTerr.setAsParent(iter, false);
        end
    end
    
    % Update handles structure
    guidata(hObject, handles);
    


function refineButton_Callback(hObject, eventdata, handles)
    pIDList = handles.evoTerr.parentIDList();
    hParentMod = cell(1,handles.evoTerr.parentCount);
    parentMod_hData = cell(1,handles.evoTerr.parentCount);
    
    for iterI=1:handles.evoTerr.parentCount
        hParentMod{iterI} = ParentModGUI(pIDList(iterI));
        parentMod_hData{iterI} = guidata(hParentMod{iterI});
        parentMod_hData{iterI}.evoGUI_handle = handles;
        guidata(hParentMod{iterI}, parentMod_hData{iterI});
    end
