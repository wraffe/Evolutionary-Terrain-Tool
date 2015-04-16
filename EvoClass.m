classdef EvoClass < handle
    %EVOCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (SetAccess = private)
        popSize;
        crossoverRate;
        globalMutationRate;
        parentLimit;
        
        terrRes;        
        sqrNumPatches; %Number of patches in a single row or column in the final rendered terrains
        overlapPerc; %Percent of patch to use in an overlap band
        stitchMethod;
        
        parentVals;
        parentCount;
        sampleFiles;
        heightMaps;
        patchMaps;
        parentPatchFitness;
        patchDatabase;
    end
    
    properties (Dependent = true, SetAccess = private)
       basePatchRes; %Resolution derived by (Terrain resolution / Number of patches in a row)
       fullPatchRes; %Resolution of each patch including the overlap bands on right and bottom
       overlapSize;  %Result of overlap percentage applied to basePatchRes
       parentIDList; %List of parents with a true value in parentVals array
       relativeMutationRate; %As the number of selected patches increases, so does the mutation rate 
    end
    
    
    methods  
        %% Constructors and Initialisers
        function EC = EvoClass(sampleTerrs, terrainResolution, squareNumPatches, patchOverlapPerc, stitchMeth, populationSize, crossRate, mutateRate)
            % Default values if constructor is called with no arguements
            EC.popSize = 8;
            EC.parentLimit = 2;
            EC.crossoverRate = 0.7;
            EC.globalMutationRate = 0.1;
            
            EC.sqrNumPatches = 4;
            EC.overlapPerc = 0.5;           
            EC.terrRes = 513;
            EC.stitchMethod = 'linear';            
            
            % Construction if arguements are provided
            if nargin > 0
                EC.terrRes = terrainResolution;
                EC.sqrNumPatches = squareNumPatches;
                EC.overlapPerc = patchOverlapPerc;
                EC.stitchMethod = stitchMeth;
                EC.popSize = populationSize;
                EC.crossoverRate = crossRate;
                EC.globalMutationRate = mutateRate;
                
                if ~isempty(sampleTerrs)                    
                    for iter = 1:size(sampleTerrs,2)
                        EC.sampleFiles{iter} = sampleTerrs{iter}; 
                    end                     
                else
                    for iter = 1:EC.popSize
                       EC.sampleFiles{iter} = strcat('heightMaps\\heightmap', int2str(iter), '.txt'); 
                    end
                end
            else
                for iter = 1:EC.popSize
                    EC.sampleFiles{iter} = strcat('heightMaps\\heightmap', int2str(iter), '.txt'); 
                end
            end
            
            EC.heightMaps = cell(1,EC.popSize);
            EC.patchMaps = cell(1,EC.popSize); 
            EC.parentPatchFitness = cell(1, EC.popSize);
            
            EC.parentVals = zeros(1,EC.popSize);
            EC.parentCount = 0;
            
            for iter = 1:EC.popSize
                EC.heightMaps{iter} = zeros(EC.terrRes, EC.terrRes); 
                EC.patchMaps{iter} = zeros(EC.sqrNumPatches, EC.sqrNumPatches);
                EC.parentPatchFitness{iter} = false(EC.sqrNumPatches, EC.sqrNumPatches);
            end
        end
        
        
        function initEvo(EC)
            %Change rand seed to be based off clock, makes sure rand seed
            %is different in every execution
            RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

            % Just a filler for now to load some terrains
            patchCount = 1;
            for iter = 1:size(EC.sampleFiles,2)
                inFileName = EC.sampleFiles{iter};
                disp(inFileName);
                tempHeightMap = readHeightMatrix(inFileName);
                if size(tempHeightMap,1)~=EC.terrRes || size(tempHeightMap,2)~=EC.terrRes
                   error('Height-maps must be square and have the same resolution as the parameter in EvoClass'); 
                end                
                EC.heightMaps{iter} = tempHeightMap; 
                
                % Extract patches from the given terrain and append them 
                % to the patch database
                tempPatchList = EC.extractPatches(EC.heightMaps{iter});
                tempPatchMap = zeros(EC.sqrNumPatches);
                
                for row = 1:EC.sqrNumPatches
                    for col = 1:EC.sqrNumPatches
                        tempPatchMap(row,col) = patchCount;
                        patchCount = patchCount+1;
                    end
                end
                EC.patchDatabase = [EC.patchDatabase tempPatchList];
                EC.patchMaps{iter} = tempPatchMap;
            end    

            EC.buildAllTerr();
        end
        
        
        %% Evolution methods
        function evolve(EC)            
            if EC.parentCount == 0
                % If no parents are selected from the current population then
                % just randomize the next generation
                for iterI = 1:EC.popSize
                    EC.patchMaps{iterI} = EC.randomTerrGen(); 
                end                
            elseif EC.parentCount == 1
                % If there is only one parent selected, then just do
                % mutation on that parent
                newGeneration = cell(1,EC.popSize);
                
                % Copy parent to next generation
                pIDList = EC.parentIDList();
                newGeneration{1} = EC.patchMaps{pIDList(1,1)};
                
                % Mutate all other offspring
                for iterI = 2:EC.popSize
                    newGeneration{iterI} = EC.mutationGen(newGeneration{1}, EC.parentPatchFitness{pIDList(1,1)});
                end 
                
                EC.patchMaps = newGeneration;
            else
                % If there is more than one parent, then do both crossover 
                % and muation
                newGeneration = cell(1,EC.popSize);
                pIDList = EC.parentIDList();                
                
                for iterI = 1:EC.popSize
                    if iterI <= EC.parentCount
                        % Copy the parents to the first few slots of the next
                        % generation
                        newGeneration{iterI} = EC.patchMaps{pIDList(1,iterI)};
                    else
                        % Conduct crossover between parents to generate the rest of
                        % the patch-maps
                        [offspringPM offspringPatchFitness]= EC.crossoverGen(pIDList);
                        offspringPM = EC.mutationGen(offspringPM, offspringPatchFitness);
                        newGeneration{iterI} = offspringPM;
                    end
                end
                
                EC.patchMaps = newGeneration;
            end           
            
            % Clear all of the per patch fitness values in preparation for
            % the next round of selection
            for iterI = 1:EC.popSize
               EC.parentPatchFitness{iterI} = false(EC.sqrNumPatches, EC.sqrNumPatches); 
            end
            
            EC.buildAllTerr();
        end
        
        function guidedEvolution(EC, targetHMFile)
            % Load target height-map
            targetHM = readHeightMatrix(targetHMFile);
            if size(targetHM,1)~=EC.terrRes || size(targetHM,2)~=EC.terrRes
                error('Target height-map must be square and have the same resolution as the parameter in EvoClass'); 
            end 
            
            % Loop until fitness reached or other generation condition reached
            for generation = 1:12
                disp(generation);
                fitnesses = zeros(EC.popSize,2);
                
                % Get fitness of all population members
                for iterI = 1:EC.popSize
                    fitnesses(iterI,1) = iterI;
                    candidateHM = EC.heightMaps{iterI};
                    
                    % Total pixel difference
                    totalDist = 0;
                    for row = 1:EC.terrRes
                        for col = 1:EC.terrRes
                            dist = abs(candidateHM(row,col) - targetHM(row,col));
                            totalDist = totalDist + dist;
                        end
                    end
                    
                    % Average pixel difference
                    if totalDist ~= 0
                        fitnesses(iterI,2) = totalDist / (EC.terrRes*EC.terrRes);
                    else
                        fitnesses(iterI,2) = 0;
                    end
                end
                    
                % Select parents
                sortedFitnesses = sortrows(fitnesses,2);
                for iterJ = 1:EC.parentLimit
                    EC.setAsParent(sortedFitnesses(iterJ,1), true);
                end
                
                % Evolve new generation
                EC.evolve();
                
                % Clear parents
                for iterJ = 1:EC.parentLimit
                    EC.setAsParent(sortedFitnesses(iterJ,1), false);
                end
                
                % Render highest parent
                if mod(generation,3) == 0
                    disp('drawing!');
                    bestFitnessHM = EC.heightMaps{sortedFitnesses(1,1)};
                    surf(bestFitnessHM, 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
                    camlight('right');
                    camproj('perspective');
                    camzoom(1.8);
                    daspect([125 125 1])
                    axis('off');
                    drawnow;
                end
            end 
            
            % Draw the target height-map
            figure;
            surf(targetHM, 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
            camlight('right');
            camproj('perspective');
            camzoom(1.8);
            daspect([125 125 1])
            axis('off');
            drawnow;
        end
        
        %% Patch map altering methods
        % Generate a new patch-map by using crossover on the selected
        % parents. pIDList should be a list of parent identifier numbers,
        % each pointing to a population member to be a parent
        function [newPatchMap newPatchFitnesses] = crossoverGen(EC, pIDList)
            % Only do crossover if there are two selected parents
            if EC.parentCount < 2
                error('There must be at least two parents selected to do crossover.');
            end            
            
            % Randomly choose one of the parents to be the base layout
            % for the new patch-map
            baseTerr = 1;%randi(EC.parentCount);
            newPatchMap = EC.patchMaps{pIDList(baseTerr)};
            basePatchFitness = EC.parentPatchFitness{pIDList(baseTerr)};
            
            % Crossover only has two parents so randomly choose another
            % parent to be a mating partner
            validPartner = false;
            partnerPatchMap = zeros(EC.sqrNumPatches, EC.sqrNumPatches);
            partnerPatchFitness = false(EC.sqrNumPatches, EC.sqrNumPatches);
            
            while ~validPartner
                partnerTerr = randi(EC.parentCount);
                % Make sure the base terrain isn't mating with itself
                if partnerTerr ~= baseTerr
                    validPartner = true;
                    partnerPatchMap = EC.patchMaps{pIDList(partnerTerr)};
                    partnerPatchFitness = EC.parentPatchFitness{pIDList(partnerTerr)};
                end
            end
            
            % Create a crossover probability matrix. Each entry corresponds to a patch            
            cProbabilities = rand(EC.sqrNumPatches);
            newPatchFitnesses = false(EC.sqrNumPatches, EC.sqrNumPatches);
            
            % Go through the patches and either keep the current patch or 
            % take a patch in the same position from another parent
            for row = 1:EC.sqrNumPatches
                for col = 1:EC.sqrNumPatches
                    % Check if this patch for either parent has been
                    % favored over the other
                    if basePatchFitness(row,col) == true && partnerPatchFitness(row,col) == false
                        % Dont change patch map, just make sure mutation
                        % knows not to change this patch
                        newPatchFitnesses(row,col) = true;
                    elseif basePatchFitness(row,col) == false && partnerPatchFitness(row,col) == true
                        % Switch the patch to the partners because it is
                        % favored
                        newPatchMap(row,col) = partnerPatchMap(row,col);
                        newPatchFitnesses(row,col) = true;
                    else
                        if basePatchFitness(row,col) == true && partnerPatchFitness(row,col) == true
                            newPatchFitnesses(row,col) = true;
                        end
                        
                        % If the probability value for a patch is less than the
                        % crossover rate then that patch shoud be switched with the
                        % equivilant patch from another parent
                        if cProbabilities(row,col) <= EC.crossoverRate
                            newPatchMap(row,col) = partnerPatchMap(row,col);
                        end
                    end
                end
            end            
        end
        
        
        % Mutate the patch map by randomly choosing patches to switch for
        % other randomly chosen patches from the database
        function [newPatchMap] = mutationGen(EC, pPatchMap, pPatchFitness)    
            % Create a mutation probability matrix. Each entry corresponds to a patch 
            mProbabilities = rand(EC.sqrNumPatches);
            
            % Set the base parent patch-map
            newPatchMap = pPatchMap;
            
            for row = 1:EC.sqrNumPatches
                for col = 1:EC.sqrNumPatches
                    if pPatchFitness(row,col) == false
                        % If the this patch is to mutate, then replace it with
                        % a randomly chosen patch from the database. Otherwise,
                        % keep this patch.
                        if mProbabilities(row,col) < EC.relativeMutationRate % Use mutation based on how many patches are being affected
                            randPatchID = randi(size(EC.patchDatabase));                       
                            newPatchMap(row,col) = randPatchID;
                        end                    
                    else
                        %disp('Patch Kept');
                    end                        
                end
            end
        end
        
        function [newPatchMap] = randomTerrGen(EC)            
            newPatchMap = randi(size(EC.patchDatabase), EC.sqrNumPatches); %make one extra for end 
            %newPatchMap = [ 17 22 33 44 68; 21 5 32 9 21; 15 6 52 36 73; 20 3 38 61 16; 1 8 35 78 43];
            %newPatchMap = [11 22 22 11 11; 11 11 11 11 11; 11 11 11 11 11; 11 11 11 11 11; 11 11 11 11 11];
        end % randomTerrGen function
        
        %% Utility terrain methods
        
        % Build a height-map from the given patch-map
        function [newTerr] = buildTerr(EC, patchMap)
            newTerr = zeros(EC.terrRes,EC.terrRes); 
                        
            fullPatchRes = EC.fullPatchRes;
            basePatchRes = EC.basePatchRes;
            
            % Build terrain starting from top left corner and moving
            % through it patch by patch, row by row
            vertPos = 1; % vertPos is bound to row loop
            for row = 0:EC.sqrNumPatches-1
                horPos = 1; % horPos is bound to column loop 
                rowMatrix = zeros(fullPatchRes, EC.terrRes);
                
                % Stitch a full row together
                for col = 0:EC.sqrNumPatches-1
                    currPatch = EC.patchDatabase{patchMap(row+1,col+1)};                    
                    
                    if col == 0
                        % Just place the first patch, no sowing or 
                        % overlaping is needed
                        rowMatrix(:, 1:fullPatchRes) = currPatch(:,:);
                        horPos = horPos + fullPatchRes;
                    else % Otherwise, stitch the new patch to the rest of the row matrix                        
                        % Get a patch worth of data to the imediate left 
                        % of where the new patch is going row matrix.
                        leftPatch = rowMatrix(:, horPos-fullPatchRes:horPos-1);
                        
                        % Stitch the new patch to the right of the existing
                        % patch
                        if strcmp(EC.stitchMethod, 'linear')
                            stitchedPatches = EC.LinearInterpStitch(leftPatch, currPatch, 'right');
                        elseif strcmp(EC.stitchMethod, 'spline')
                            stitchedPatches = EC.CubicSplineStitch(leftPatch, currPatch, 'right');
                        else
                            error('Invalid stitching method parameter supplied');
                        end                        
                        
                        % If the whole stitched result wont fit back into
                        % the row matrix, just copy as much as possible.
                        % Else, just copy the whole stitched result back in
                        if horPos+basePatchRes-1 > size(rowMatrix,2)
                            pixelsLeft = size(rowMatrix,2)-(horPos-fullPatchRes)+1;
                            rowMatrix(:, horPos-fullPatchRes:size(rowMatrix,2)) = stitchedPatches(:, 1:pixelsLeft);        
                        else
                            rowMatrix(:, horPos-fullPatchRes:horPos+basePatchRes-1) = stitchedPatches(:,:);                            
                        end
                        
                        % Position isn't incremented by full patch res
                        % because overlapping reduces the size of the new patch
                        horPos = horPos + basePatchRes;                        
                    end % row stitching conditions                   
                end % col loop                 
                
                % Put the newly created row into the final terrain
                if row == 0
                    % If this is the first row, simply copy it into the
                    % terrain
                    newTerr(vertPos:vertPos+fullPatchRes-1, :) = rowMatrix(:,:); 
                    vertPos =  vertPos + fullPatchRes;
                else
                    % Else, stitch the row onto the previous row that is
                    % already in the terrain
                    topPatch = newTerr(vertPos-fullPatchRes:vertPos-1,:);
                    if strcmp(EC.stitchMethod, 'linear')
                        stitchedRows = EC.LinearInterpStitch(topPatch, rowMatrix, 'below');
                    elseif strcmp(EC.stitchMethod, 'spline')
                        stitchedRows = EC.CubicSplineStitch(topPatch, rowMatrix, 'below');
                    else
                        error('Invalid stitching method parameter supplied');
                    end 
                    
                    % If its the last row, dont copy all of the content,
                    % only what can fit
                    if vertPos+basePatchRes-1 > size(newTerr,1)
                        pixelsLeft = size(newTerr,1)-(vertPos-fullPatchRes)+1;
                        newTerr(vertPos-fullPatchRes:size(newTerr,1),:) = stitchedRows(1:pixelsLeft,:);
                    else
                        newTerr(vertPos-fullPatchRes:vertPos+basePatchRes-1,:) = stitchedRows(:,:);
                    end
                    
                    vertPos = vertPos + basePatchRes;
                end                
            end % row loop 
        end % End buildTerr function
        
        
        % Batch build all of the terrains in the population from their
        % corresponding patch-map
        function buildAllTerr(EC)
           for iterI = 1:EC.popSize
              EC.heightMaps{iterI} = EC.buildTerr(EC.patchMaps{iterI}); 
           end
        end
        
        
        function [stitchedPatches] = LinearInterpStitch(EC, fixedPatch, newPatch, stitchDir)    
            if (size(fixedPatch,1) ~= size(newPatch,1) || size(fixedPatch,2) ~= size(newPatch,2))
                error('Error in stitching: Patches must have the same dimensions');
            end

            %patchRes = size(newPatch,1);

            % Defining these seperately allows for non-square patches to be
            % supplied
            patchHeight = size(newPatch,1);
            patchLength = size(newPatch,2);
            
            overlapSize = EC.overlapSize;
            %rearOverlap = EC.fullPatchRes - overlapSize + 1;

            % stitchDir is the position of newPatch relative to fixedPatch.
            % Thus, 'below' is 'newPatch is below fixedPatch' and 'right' is
            % 'newPatch is to the right of fixedPatch'            
            if strcmpi(stitchDir,'right') || strcmpi(stitchDir,'left')
                % Depending on which side the new patch is to be stitched onto the
                % existing patch, just reverse the order that patches are given to 
                % the following stitching algorithm. I.e. the below algorithm will
                % always put patch2 to the right of patch1
                if strcmpi(stitchDir,'right')
                    patch1 = fixedPatch;
                    patch2 = newPatch;
                elseif strcmpi(stitchDir,'left')
                    patch1 = newPatch;
                    patch2 = fixedPatch;
                end

                %overlapSize = EC.overlapPerc * patchLength;
                rearOverlap = patchLength - overlapSize + 1;
                
                overlap1 = patch1(:, rearOverlap:patchLength); %Copy all rows of the rear columns
                overlap2 = patch2(:, 1:overlapSize); %Copy all rows of the leading columns

                overlapsMerged = zeros(patchHeight,overlapSize);
                for iterJ=1:overlapSize
                    percentage = iterJ / overlapSize;

                    for iterI=1:patchHeight
                        left = overlap1(iterI,iterJ);
                        right = overlap2(iterI,iterJ);

                        dist = right - left;
                        linear = dist * percentage;

                        overlapsMerged(iterI,iterJ) = left + linear;
                    end
                end

                stitchedSize = patchLength*2-overlapSize;
                stitchedPatches = zeros(patchHeight,stitchedSize);
                stitchedPatches(:,1:rearOverlap-1) = patch1(:, 1:rearOverlap-1);
                stitchedPatches(:,rearOverlap:patchLength) = overlapsMerged(:,1:overlapSize);
                stitchedPatches(:,patchLength+1:stitchedSize) = patch2(:,overlapSize+1:patchLength); 

            elseif strcmpi(stitchDir,'below') || strcmpi(stitchDir,'above') 
                % patch2 is always below patch1. Again, decide where incoming
                % patches are assigned
                if strcmpi(stitchDir,'below')
                    patch1 = fixedPatch;
                    patch2 = newPatch;
                elseif strcmpi(stitchDir,'above')
                    patch1 = newPatch;
                    patch2 = fixedPatch;
                end

                %overlapSize = EC.overlapPerc * patchHeight;
                rearOverlap = patchHeight - overlapSize + 1;
                
                overlap1 = patch1(rearOverlap:patchHeight, :); %Copy all rows of the rear columns
                overlap2 = patch2(1:overlapSize, :); %Copy all rows of the leading columns

                overlapsMerged = zeros(overlapSize,patchLength);
                for iterJ=1:overlapSize
                    percentage = iterJ / overlapSize;

                    for iterI=1:patchLength
                        top = overlap1(iterJ,iterI);
                        bottom = overlap2(iterJ,iterI);

                        dist = bottom - top;
                        linear = dist * percentage;

                        overlapsMerged(iterJ,iterI) = top + linear;
                    end
                end

                stitchedSize = patchHeight*2-overlapSize;
                stitchedPatches = zeros(stitchedSize,patchLength);
                stitchedPatches(1:rearOverlap-1,:) = patch1(1:rearOverlap-1,:);
                stitchedPatches(rearOverlap:patchHeight,:) = overlapsMerged(1:overlapSize,:);
                stitchedPatches(patchHeight+1:stitchedSize,:) = patch2(overlapSize+1:patchHeight,:);   
            else        
                error('The patch stitching direction provided is not valid');
            end            
        end
        
        function [stitchedPatches] = CubicSplineStitch(EC, fixedPatch, newPatch, stitchDir)
            if (size(fixedPatch,1) ~= size(newPatch,1) || size(fixedPatch,2) ~= size(newPatch,2))
                error('Error in stitching: Patches must have the same dimensions');
            end

            % Defining these seperately allows for non-square patches to be
            % supplied
            patchHeight = size(newPatch,1);
            patchLength = size(newPatch,2);
            
            overlapSize = EC.overlapSize;
            
            % Create a spline curve to get the percentage values to
            % interpolate upon. Avoids linear interpolation.            
            x = [0 0.1 0.25 0.5 0.75 0.9 1]; % x positions of control points
            y = [0 0.02 0.15 0.5 0.85 0.98 1]; % y positions of control points
            cSpline = spline(x,y);
            %disp(ppval(cs,1));
            %xx = linspace(1,4,100);
            %plot(x,y,'o',xx,ppval(cs,xx),'-');

            % stitchDir is the position of newPatch relative to fixedPatch.
            % Thus, 'below' is 'newPatch is below fixedPatch' and 'right' is
            % 'newPatch is to the right of fixedPatch'            
            if strcmpi(stitchDir,'right') || strcmpi(stitchDir,'left')
                % Depending on which side the new patch is to be stitched onto the
                % existing patch, just reverse the order that patches are given to 
                % the following stitching algorithm. I.e. the below algorithm will
                % always put patch2 to the right of patch1
                if strcmpi(stitchDir,'right')
                    patch1 = fixedPatch;
                    patch2 = newPatch;
                elseif strcmpi(stitchDir,'left')
                    patch1 = newPatch;
                    patch2 = fixedPatch;
                end

                rearOverlap = patchLength - overlapSize + 1;
                
                overlap1 = patch1(:, rearOverlap:patchLength); %Copy all rows of the rear columns
                overlap2 = patch2(:, 1:overlapSize); %Copy all rows of the leading columns

                overlapsMerged = zeros(patchHeight,overlapSize);
                for iterJ=1:overlapSize
                    % How far along the spline curve to get our percentage
                    % value
                    interpDist = iterJ / overlapSize; 
                    percentage = ppval(cSpline,interpDist);

                    for iterI=1:patchHeight
                        left = overlap1(iterI,iterJ);
                        right = overlap2(iterI,iterJ);

                        dist = right - left;
                        linear = dist * percentage;

                        overlapsMerged(iterI,iterJ) = left + linear;
                    end
                end

                stitchedSize = patchLength*2-overlapSize;
                stitchedPatches = zeros(patchHeight,stitchedSize);
                stitchedPatches(:,1:rearOverlap-1) = patch1(:, 1:rearOverlap-1);
                stitchedPatches(:,rearOverlap:patchLength) = overlapsMerged(:,1:overlapSize);
                stitchedPatches(:,patchLength+1:stitchedSize) = patch2(:,overlapSize+1:patchLength); 

            elseif strcmpi(stitchDir,'below') || strcmpi(stitchDir,'above') 
                % patch2 is always below patch1. Again, decide where incoming
                % patches are assigned
                if strcmpi(stitchDir,'below')
                    patch1 = fixedPatch;
                    patch2 = newPatch;
                elseif strcmpi(stitchDir,'above')
                    patch1 = newPatch;
                    patch2 = fixedPatch;
                end

                rearOverlap = patchHeight - overlapSize + 1;
                
                overlap1 = patch1(rearOverlap:patchHeight, :); %Copy all rows of the rear columns
                overlap2 = patch2(1:overlapSize, :); %Copy all rows of the leading columns

                overlapsMerged = zeros(overlapSize,patchLength);
                for iterJ=1:overlapSize
                    % How far along the spline curve to get our percentage
                    % value
                    interpDist = iterJ / overlapSize; 
                    percentage = ppval(cSpline,interpDist);

                    for iterI=1:patchLength
                        top = overlap1(iterJ,iterI);
                        bottom = overlap2(iterJ,iterI);

                        dist = bottom - top;
                        linear = dist * percentage;

                        overlapsMerged(iterJ,iterI) = top + linear;
                    end
                end

                stitchedSize = patchHeight*2-overlapSize;
                stitchedPatches = zeros(stitchedSize,patchLength);
                stitchedPatches(1:rearOverlap-1,:) = patch1(1:rearOverlap-1,:);
                stitchedPatches(rearOverlap:patchHeight,:) = overlapsMerged(1:overlapSize,:);
                stitchedPatches(patchHeight+1:stitchedSize,:) = patch2(overlapSize+1:patchHeight,:);   
            else        
                error('The patch stitching direction provided is not valid');
            end       
        end
        
        function [patchArray] = extractPatches(EC, fullHM)
            fullPatchRes = EC.fullPatchRes;
            basePatchRes = EC.basePatchRes;
            
            fullHMRes = size(fullHM,1);     
    
            nPatchesInRow = floor(fullHMRes / basePatchRes);
    
            patchArray = cell(1,(nPatchesInRow*nPatchesInRow));
    
            numPatches = 0;
            rowPt = 1;
    
            for rowOffset = 0:nPatchesInRow-1
                nextRowPt = rowPt+fullPatchRes-1;

                colPt = 1;           
                for colOffset = 0:nPatchesInRow-1
                    nextColPt = colPt+fullPatchRes-1;
                    numPatches = numPatches+1;
                    
                    if nextColPt > fullHMRes && nextRowPt < fullHMRes % If last column
                        % If this is the last patch in the row, then
                        % unfold it to make it a full sized patch
                        hPixelsLeft = fullHMRes - colPt;
                        subPatch = fullHM(rowPt:nextRowPt, colPt:colPt+hPixelsLeft);
                        newPatch = EC.unfoldPatch(subPatch);
                    elseif nextColPt < fullHMRes && nextRowPt > fullHMRes % If last row
                        vPixelsLeft = fullHMRes - rowPt;
                        subPatch = fullHM(rowPt:rowPt+vPixelsLeft, colPt:nextColPt);
                        newPatch = EC.unfoldPatch(subPatch);
                    elseif nextColPt > fullHMRes && nextRowPt > fullHMRes % If last patch                        
                        hPixelsLeft = fullHMRes - colPt;
                        vPixelsLeft = fullHMRes - rowPt;
                        subPatch = fullHM(rowPt:rowPt+vPixelsLeft, colPt:colPt+hPixelsLeft);
                        newPatch = EC.unfoldPatch(subPatch);
                    else % If a normal center patch
                        newPatch = fullHM(rowPt:nextRowPt, colPt:nextColPt);
                    end                    
                    
                    patchArray{numPatches} =  newPatch;                    
                        
                    colPt = colPt+basePatchRes;
                end
                
                rowPt = rowPt+basePatchRes;
            end    
            
            % Delete cells so that the cell array has no blank slots
            %numPatches = numPatches+1;
            %while (numPatches < size(patchArray,2))
                %patchArray(numPatches) = [];
            %end
        end  
        
        
        % Given a small patch, expand it out to a full sized patch by
        % mirroring it constantly
        function [finalPatch] = unfoldPatch(EC, subPatch)
            fullPatchRes = EC.fullPatchRes;

            finalPatch = zeros(fullPatchRes, fullPatchRes);            
            
            vPos = size(subPatch,1);
            hPos = size(subPatch,2); 
            
            finalPatch(1:vPos,1:hPos) = subPatch(:,:);
           
            % Flip horizontally first, and then verticaly            
            hFlipSubPatch = subPatch; 
            while hPos < fullPatchRes                
                endHPos = hPos + size(subPatch,2) - 1;
                hFlipSubPatch = fliplr(hFlipSubPatch);                
                
                if endHPos > fullPatchRes  
                    % If cant coppy full amount, only coppy as much as will
                    % fit
                    pixelsLeft = fullPatchRes - hPos;
                    finalPatch(1:vPos, hPos:hPos+pixelsLeft) = hFlipSubPatch(:,1:pixelsLeft+1);
                else
                    % Else, copy the whole thing
                    finalPatch(1:vPos, hPos:endHPos) = hFlipSubPatch(:,:);
                end
                
                hPos = hPos + size(subPatch,2);
            end
            
            % Take the full row of the final patch to keep mirroring
            vFlipSubPatch = finalPatch(1:vPos,:);
            while vPos < fullPatchRes               
               endVPos = vPos + size(subPatch,1) - 1;
               vFlipSubPatch = flipud(vFlipSubPatch);
               
               if endVPos > fullPatchRes
                    pixelsLeft = fullPatchRes - vPos;
                    finalPatch(vPos:vPos+pixelsLeft, :) = vFlipSubPatch(1:pixelsLeft+1,:);
               else
                    finalPatch(vPos:endVPos, :) = vFlipSubPatch(:,:); 
               end
               vPos = vPos + size(subPatch,1);
            end
        end     
        
        
        function patchShow(EC, terrainNum)
            patchMap = EC.patchMaps{terrainNum};
            heightMap = EC.heightMaps{terrainNum};
            
            figure;
            surf(heightMap, 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
            view(0,35);
            camlight('right');
            camproj('perspective');
            camzoom(1.8);
            daspect([125 125 1])
            axis('off');
            cRatio = caxis; 
            
            for row = 1:EC.sqrNumPatches
                for col = 1:EC.sqrNumPatches                     
                    patchID = patchMap(row,col);                    
                    patchData = EC.patchDatabase{patchID};
                    figure;
                    surf(patchData, 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
                    caxis(cRatio)
                    view(0,35);
                    camlight('right');
                    camproj('orthographic');%perspective');
                    camzoom(0.5);
                    daspect([125 125 1])                    
                    axis('off');
                end
            end           
        end
        
        
        function [markedHM] = markedPatchBoarders(EC, terrainNum)
            markedHM = EC.heightMaps{terrainNum};
            basePatchRes= EC.basePatchRes();
            boarderColor = -1;
            
            for row = 1:EC.terrRes
                for col = 1:EC.terrRes
                    rowRemainder = mod(row, basePatchRes);
                    colRemainder = mod(col, basePatchRes);                   
                    
                    if colRemainder == 0 || rowRemainder == 0 || row == 5 || col == 5                        
                        markedHM(row,col) = boarderColor; 
                    end
                end
            end
        end
        
        %% Set methods
        function [success] = setAsParent(EC, terrID, bool)           
            if terrID < 1 || terrID > EC.popSize
                error('Terrain ID is not valid. Terrain ID must be between 1 and the population size parameter');               
            end            
            
            if bool == true
                if EC.parentCount >= EC.parentLimit
                    success = false;
                else
                    EC.parentVals(terrID) = true;
                    EC.parentCount = EC.parentCount + 1;
                    success = true;
                end
            elseif bool == false
                EC.parentVals(terrID) = false;
                if EC.parentCount > 0
                    EC.parentCount = EC.parentCount - 1;
                end
                success = true;
            else
                error('Must provide a true or false for bool in setAsParent()');
            end            
        end
        
        
        function setParentPatchFitness(EC, hmID, row, col, bool)
            if hmID < 1 || hmID > EC.popSize
                error('Terrain ID is not valid. Terrain ID must be between 1 and the population size parameter');               
            end    
            if row > EC.sqrNumPatches ||  row < 0 || col > EC.sqrNumPatches || col < 0
                error('Row and column identifiers must be between 0 and the sqrNumPatches parameter');
            end
            
            tempMap = EC.parentPatchFitness{hmID};
            tempMap(row,col) = bool;
            EC.parentPatchFitness{hmID} = tempMap;
        end     
        
        %% Get methods
        function [matrix] = getHeightMap(EC, popID)
            matrix = EC.heightMaps{popID};
        end       
        
        function [matrix] = getPatchMap(EC, popID)
            matrix =  EC.patchMaps{popID};
        end
        
        function [matrix] = getParentPatchFitness(EC, popID)
           matrix = EC.parentPatchFitness{popID}; 
        end
        
        function [patch] = getPatch(EC, patchID)
           patch = EC.patchDatabase{patchID};
        end           
        
        function [patchResolution] = get.basePatchRes(EC)
            patchResolution = floor(EC.terrRes / EC.sqrNumPatches);       
        end
        
        function [patchNOverlap] = get.fullPatchRes(EC)
            patchRes = EC.basePatchRes;
            patchNOverlap = patchRes + EC.overlapSize;
        end        
        
        function [size] = get.overlapSize(EC)
            % overlapPerc is percent of patch to overlap with an adjacent patch.
            % Thus, if overlapPerc = 0.25, then 25% of the patch will overlap with
            % patches on either side of this patch, thus 50% of columns in patch
            % will be used for overlapping in a center patch. overlapSize is number 
            % of pixels in each overlap area, rounded down.
            size = floor(EC.basePatchRes*EC.overlapPerc);
        end 
        
        function [pList] = get.parentIDList(EC)
            pCount = 0;
            pList = zeros(1, EC.parentCount);
            for iterI = 1:EC.popSize
                if EC.parentVals(iterI) == true
                    pCount = pCount+1;
                    pList(pCount) = iterI;
                end
            end
        end
        
        function [finalMutationRate] = get.relativeMutationRate(EC)
           % As the number of patches selected during the refining parent 
           % stage increases, so does the mutation rate. This means that if
           % there are fewer patches that are allowed to mutate then they
           % are for more likely to do so. This prevents entire generations
           % of identical terrains being generate when the user only wants
           % a few patches to change
           
           % Get average number of selected parent patches from all parents
           selectedPatchCount = 0;
           pIDList = EC.parentIDList;
           for iterI = 1:EC.parentCount
               patchFitnessMap = EC.getParentPatchFitness(pIDList(iterI));
               for row = 1:EC.sqrNumPatches
                  for col = 1:EC.sqrNumPatches
                      if patchFitnessMap(row,col) == true
                          selectedPatchCount = selectedPatchCount + 1;
                      end
                  end
              end
           end
           
           avgPatchCount = floor(selectedPatchCount/EC.parentCount);
           
           percentPatchesSelected = avgPatchCount / (EC.sqrNumPatches*EC.sqrNumPatches);
           
           mutationRange = 1-EC.globalMutationRate;
           mutationPercentIncrease = percentPatchesSelected * mutationRange;
           finalMutationRate = EC.globalMutationRate + mutationPercentIncrease;
           
        end
    end %End public methods  
end %End class












