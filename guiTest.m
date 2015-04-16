clc;
clear;


sampleFiles = cell(1,8);
for iter = 1:size(sampleFiles,2)
    sampleFiles{iter} = strcat('heightMaps\\heightmap', int2str(iter), '.txt');
end

%EvoGUI;
%ParentModGUI;

evo = EvoClass(sampleFiles,513,2,0.5,'linear',4,0.5,0.5); 
evo.initEvo();
%evo.patchShow(1);
%hm = evo.getHeightMap(8);
%patch = hm(201:250, 1:50);
%expandedPatch = evo.unfoldPatch(patch);
%imshow(expandedPatch);
%evo.evolve();
%evo.guidedEvolution('heightMaps\\targetHM.txt');
terrain = evo.getHeightMap(8);
%patch1 = evo.getPatch(18);
%patch2 = evo.getPatch(14);
%stitched = linearPatchStitch(patch1,patch2,'left',0.5); 
%overlapSz = evo.overlapSize;
%stitched = SplinePatchStitch(overlapSz,patch1,patch2,'above');


%cm = [0.839 0.714 0.605]; %pink
%colormap(cm);
surf(terrain, 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
caxis([-0.01 1]);
%view(0,40);
camlight('right');
%camlight('right');
camproj('perspective');
camzoom(1.8);
daspect([125 125 1])
axis('off');



