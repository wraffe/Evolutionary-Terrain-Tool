function varargout = ParentModGUI(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ParentModGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ParentModGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


function ParentModGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ParentModGUI
handles.output = hObject;

set(handles.topDownAxis, 'Visible', 'off');
set(handles.perspectiveAxis, 'Visible', 'off');

handles.hmID = varargin{1};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ParentModGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function varargout = ParentModGUI_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;


function closeButton_Callback(hObject, eventdata, handles)
    delete(handles.figure1);

function loadButton_Callback(hObject, eventdata, handles)
    evoGUIHandle = handles.evoGUI_handle;
    evoClass = evoGUIHandle.evoTerr;
    parentHM = evoClass.getHeightMap(handles.hmID);
    % Get height-map with patch grids marked on it
    markedParentHM = evoClass.markedPatchBoarders(handles.hmID);
    
    % Draw terrain
    set(handles.topDownAxis, 'Visible', 'on');
    set(handles.perspectiveAxis, 'Visible', 'on');

    % Draw top down view axis
    surf(handles.topDownAxis, markedParentHM, 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    caxis(handles.topDownAxis,[0 1]);
    cm = [0.839 0.714 0.605]; %pink
    colormap(handles.topDownAxis, cm);
    %view(handles.topDownAxis, 180, 90);
    view(handles.topDownAxis, 0, 90);
    light('Parent', handles.topDownAxis); %specifically create a light for each object
    camlight('right');
    camproj(handles.topDownAxis, 'perspective');
    %camzoom(handles.topDownAxis, 1.8);
    daspect(handles.topDownAxis, [125 125 1])
    axis(handles.topDownAxis, 'off');
    xlim(handles.topDownAxis, [1 evoClass.terrRes]); % Prevents extra space being taken up by a larger grid than is necessary
    ylim(handles.topDownAxis, [1 evoClass.terrRes]);
        
    % Draw perspective view axis
    surf(handles.perspectiveAxis, parentHM, 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    caxis(handles.perspectiveAxis,[0 1]);
    cm = [0.839 0.714 0.605]; %pink
    colormap(handles.perspectiveAxis, cm);
    %view(handles.perspectiveAxis, 140, 30);
    view(handles.perspectiveAxis, 0, 30);
    light('Parent', handles.perspectiveAxis); %specifically create a light for each object
    %camlight('right')
    camproj(handles.perspectiveAxis, 'perspective');
    camzoom(handles.perspectiveAxis, 1.8);
    daspect(handles.perspectiveAxis, [125 125 1])
    axis(handles.perspectiveAxis, 'off');
    
    
    % Create the checkbox objects
    patchFitnesses = evoClass.getParentPatchFitness(handles.hmID);
    % Each square in GUIDE is 20 by 20 (pixels?). [0,0] is bottom left of
    % window. topDownAxis is from [60,130] (bottom left) to [500,570] top right
    sqrNumPatches = evoClass.sqrNumPatches;
    boxSize = 15;
    % Takes into account the axis lines that are invisible but still take up space on gui
    axisOffset = 12; 
    posVect = get(handles.topDownAxis, 'Position');    
    width = posVect(3)-axisOffset; % width is 3rd value arguement
    dist = (width / sqrNumPatches);
    firstBoxPos = [(posVect(1)+(dist/2)) (posVect(2)+(dist/2))];
    
    idCount = 1;
    for col = 1:sqrNumPatches   
        for row = 1:sqrNumPatches
            rowOffset = (row-1)*dist;
            colOffset = (col-1)*dist;
            rowPos = firstBoxPos(1) + rowOffset;
            colPos = firstBoxPos(2) + colOffset;            
            id = idCount;
            idCount = idCount+1;  
            
            if patchFitnesses(col,row) == true
                uicontrol('Style','checkbox','Tag',int2str(id),'Value',1,'Callback',{@checkbox_Callback,handles},'String',int2str(id),'Position',[rowPos colPos boxSize boxSize]);
            else
                uicontrol('Style','checkbox','Tag',int2str(id),'Callback',{@checkbox_Callback,handles},'String',int2str(id),'Position',[rowPos colPos boxSize boxSize]);
            end
        end
    end 

    % Update handles structure
    guidata(hObject, handles);
    
    
    
function checkbox_Callback(hObject, eventData, handles)
    evoGUIHandle = handles.evoGUI_handle;
    evoClass = evoGUIHandle.evoTerr;
    
    % This callback is used by all checkboxes so we need to get the id of
    % the box and turn it into patch-map coordinates
    sqrNumPatches = evoClass.sqrNumPatches;
    id = sscanf(get(hObject,'Tag'), '%d');
    rowPos = ceil(id/sqrNumPatches);
    colPos = mod(id,sqrNumPatches);
    if colPos == 0
        colPos = 5; 
    end

    if get(hObject, 'Value') == get(hObject, 'Max')
        % If check box is ticked
        evoClass.setParentPatchFitness(handles.hmID, rowPos, colPos, true);
    else
        evoClass.setParentPatchFitness(handles.hmID, rowPos, colPos, false);
    end
