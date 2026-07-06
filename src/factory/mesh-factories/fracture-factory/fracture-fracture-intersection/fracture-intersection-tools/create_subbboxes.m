%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p1,p2,p3,p4] = create_subbboxes(diagonal1, diagonal2, N)
 % Calculate the minimum and maximum coordinates of the bounding box
    minX = min(diagonal1(1), diagonal2(1));
    maxX = max(diagonal1(1), diagonal2(1));
    minY = min(diagonal1(2), diagonal2(2));
    maxY = max(diagonal1(2), diagonal2(2));
    
    % Calculate the size of each smaller box
    boxWidth = (maxX - minX) / N;
    boxHeight = (maxY - minY) / N;
    
    % Initialize an empty cell array to store the subdivided boxes
    subdividedBoxes = cell(N, N);
    
    % Subdivide the bounding box into smaller boxes
    for i = 1:N
        for j = 1:N
            % Calculate the coordinates of the current smaller box
            boxMinX = minX + (i - 1) * boxWidth;
            boxMaxX = minX + i * boxWidth;
            boxMinY = minY + (j - 1) * boxHeight;
            boxMaxY = minY + j * boxHeight;
            
            % Store the coordinates of the current smaller box
            subdividedBoxes{i, j} = [boxMinX, boxMinY, ...
                                     boxMinX, boxMaxY, ...
                                     boxMaxX, boxMaxY, ...
                                     boxMaxX, boxMinY];
        end
    end
    
    % Calculate the coordinates of the divided bounding box
    auxmat = cell2mat(reshape(subdividedBoxes,[],1));
    p1 = auxmat(:,1:2);
    p2 = auxmat(:,3:4);
    p3 = auxmat(:,5:6);
    p4 = auxmat(:,7:8);
end
