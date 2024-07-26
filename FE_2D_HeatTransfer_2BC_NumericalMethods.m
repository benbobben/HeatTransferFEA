%% Final Assignement Neumerical Methods 

% Final Heat Flow Assignment 
% Goal: keep the room at a temperature of 21 degrees C on the coldest day. 
% Given:  of 10 degrees C and that the walls exposed to the outside (windows) have a temperature of -23 degrees C:
% Add 5 Heaters to the porch area

% l = the height of the room
% T = temperature
% Q = heat energy generated per unit area per unit time W atts = Joule=secn
% K = thermal conductivity W= (m deg C):
% (xi; yi) is the Dirac delta function
close all
clc 
clear

addpath("/Users/benjackson/Desktop/School/Grad -Spring/Neumerical Methods/Final/Mesh /distmesh-master/")

% q = W 
% k = W/mC =Thermal conductivity

%% Defining the parameters associated with the Model (K)
minVal      = .001;        % minimum value for rounding error --> 0 

thickness_window_in = 0.5;                      % in
thickness_window = thickness_window_in*0.0254;  % m

R_windows   = 1.38;                               %W/(m·C); %window 
k_window    = .96;                                % W / (mC) thermal Conductivity of Window
% k_window    = .96/60;                                % W / (mC) thermal Conductivity of Window
k_wall      = .17;                                % W / (mC) 
k_air       = .026;                               % W / (mC) Thermal conductivity of air

% Define the boundary temperatures
boundaryTemperature = -23;                      % deg C
houseWallTemperature = 10;                      % deg C % East Wall of porch

% room geometry % Assuming height of the room is constant throughout so we negelct it.
wallLength1_in = 45*3 + 3*2; %in
wallLength2_in = 45*4 + 3*3; %in
wallLength1 = wallLength1_in/39.37; %m
wallLength2 = wallLength2_in/39.37; %m

% Define dimensions of the rectangle
x_length = wallLength2; % Length along x-axis
y_length = wallLength1; % Length along y-axis

% Define the coordinates of the rectangle
x_rect = [0, x_length, x_length, 0, 0]; % x-coordinates of the rectangle
y_rect = [0, 0, y_length, y_length, 0]; % y-coordinates of the rectangle

% Plot the rectangle
figure;
plot(x_rect, y_rect, 'b', 'LineWidth', 2);
axis equal; % Make the aspect ratio equal
title('Cross Section of Room');
xlabel('X (m)');
ylabel('Y (m)');
grid on;
hold on 
bufferSize = 10/39.37 %10in from the wall /39.37 - > m
x_heated = [bufferSize    4.8006-bufferSize    4.8006-bufferSize         bufferSize         bufferSize]
y_heated = [ bufferSize         bufferSize    3.5814-bufferSize    3.5814-bufferSize        bufferSize]
plot(x_heated,y_heated, 'r', 'LineWidth', 2);
fontsize(25,"points");
hold off


%% Generate Mesh and plot
fac = 1;
fd  = @(p) -min(min(min(p(:,2),y_length-p(:,2)),0+p(:,1)),x_length-p(:,1));
fh  = @(p) ones(size(p,1),1);
[p,t] = distmesh( fd, fh, .05, [0,0;round(x_length/fac),round(y_length/fac)]); % 3rd Varibale is length of sides on elements
patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )
transformationTable = t; % t is triangle indices 
nodeLocationTable   = p; % p is  Grid vertex/node coordinates (N_P x 2/3)

%% Assign a K value for each element 
% if K is on an edge assign it either a wall or window k value 
x_max = max(nodeLocationTable(:,1));
y_max = max(nodeLocationTable(:,2));
x_min = min(nodeLocationTable(:,1));
y_min = min(nodeLocationTable(:,2));

% Window locations are everywhere except 3in from the coners and at the
% width of the wall in 3 locations 
% Defining seperate array for conductivity at each element
boundary_Array = zeros(length(nodeLocationTable),1);            % If all were same type boundary condition
boundary_Array_House = zeros(length(nodeLocationTable),1);      % Just Boundary Conditions for house facing wall
boundary_Array_Outside = zeros(length(nodeLocationTable),1);    % Just Boundary Conditions for outside facing wall
globalMatrix   = zeros(length(nodeLocationTable),length(nodeLocationTable)); %Global Matrix
k_node = zeros(length(nodeLocationTable),1);                    % Array of thermal constant at each node
k_Array = zeros(length(transformationTable),1); %An array of the size of elements to hold k value of each element
wallHeat_asFlux_Array = zeros(length(globalMatrix),1);

% What i am trying to do here is that if the current element is on the
% boundary then it will have a k value of the wall or the window
for e = 1:height(transformationTable)
    element = transformationTable(e,:);
    for i = 1:3
        %if on a boundary then k is wall or window
        node_i_X = nodeLocationTable(element(i),1);
        node_i_Y = nodeLocationTable(element(i),2);
        if node_i_X >= x_max -minVal|| node_i_Y >= y_max -minVal|| node_i_X <= minVal || node_i_Y <= minVal
            % It is a boundary element 

            %%%%%%% For not the boundary will be all windows
            k_Array(e) = k_window;
            break
        else
            % Else the element is air 
            k_Array(e) = k_air;
        end
    end
end

%% Define the boundary conditions -- The temperature at each boundary Node
%   Results section - arrays of boundary temperature 
%                   boundary_Array - which nodes are boundaries and
%                       the temperature of them - used for constant BC
%                   boundary_Array_Outside - just the outside BC
%                   boundary_Array_House - just the house facing wall
%   k_node - the node's thermal conductivity if its a wall or window                 
for node = 1:length(nodeLocationTable)  
    %If we are on the walls facing the outside temperature then set BC
    if nodeLocationTable(node,1) >= x_max - minVal || nodeLocationTable(node,1) <= minVal || nodeLocationTable(node,2) <= minVal
        boundary_Array(node)         = boundaryTemperature; % Setting the temerature at the boundary to be 20 degrees
        boundary_Array_Outside(node) = boundaryTemperature; % Setting the temerature at the boundary to be 20 degrees
    end 
    %if y = 0 on geometrey then the on the east wall and its value is houseWallTemperature
    if nodeLocationTable(node,2) >= y_max - minVal
        boundary_Array(node)       = houseWallTemperature;
        boundary_Array_House(node) = houseWallTemperature;
    end

    %Find if the current node is a window or a wall - We are on a boundary:
    if nodeLocationTable(node,1) >= x_max - minVal || nodeLocationTable(node,1) <= minVal
        % North or South Wall - If the y value on north wall is 3in from the
        % corners then its a wall 
        if nodeLocationTable(node,2) < 3/39.37 ||  nodeLocationTable(node,2) > y_max  - 3/39.37%m
            k_node(node) = k_wall;
        else 
            k_node(node) = k_window;
        end

    elseif nodeLocationTable(node,2) <= y_min + minVal
        % West Wall
        if (nodeLocationTable(node,1) < 3/39.37 ||  nodeLocationTable(node,1) > x_max - 3/39.37 || (nodeLocationTable(node,1) > x_max/2 - (3/39.37)/2) && (nodeLocationTable(node,1) < x_max/2 + (3/39.37)/2))
            k_node(node) = k_wall;
        else 
            k_node(node) = k_window;
        end
    % elseif nodeLocationTable(node,2) <= minVal
    %     % East Wall
    %     if (nodeLocationTable(node,2) < 3/39.37 ||  nodeLocationTable(node,2) > x_max - 3/39.37 || (nodeLocationTable(node,2) > 83/39.37) && (nodeLocationTable(node,2) < 86/39.37 )|| (nodeLocationTable(node,2) > (3+80+3+17)/39.37 && nodeLocationTable(node,2) < (3+80+3+17+3)/39.37))    
    %         k_node(node) = k_wall;
    %     else 
    %         k_node(node) = k_window;
    %     end
    end
end

%% Get Global Matrix - Solve Each node across the mesh 
for e = 1:height(transformationTable) % loop through the elements

    element = transformationTable(e,:); % The current element node locations

    %Get the area of the element
    P = [   nodeLocationTable(element(1),1) , nodeLocationTable(element(1),2), 1;
            nodeLocationTable(element(2),1) , nodeLocationTable(element(2),2), 1;
            nodeLocationTable(element(3),1) , nodeLocationTable(element(3),2), 1];

    A = abs(1/2 * det(P));

    % loop for the rows of the coefficient Matrix
    for i = 1:3 
        j = mod(i, 3) + 1;     % i
        k = mod(i + 1, 3) + 1; % i
        
        nodej = element(j);% The Node in from global
        nodek = element(k);% The Node in from global
        
        % Get the derivatives for Phi
        dNi_dx = (nodeLocationTable(nodej,2)   -  nodeLocationTable(nodek,2));
        dNi_dy = (nodeLocationTable(nodek,1)   -  nodeLocationTable(nodej,1));

        % loop for the columns of the coefficient Matrix
        for l = 1:3
            j = mod(l, 3)     + 1;      % j
            k = mod(l + 1, 3) + 1;  % j

            % Get derivates for Phi j
            dNj_dx = (nodeLocationTable(element(j),2)   -  nodeLocationTable(element(k),2));
            dNj_dy = (nodeLocationTable(element(k),1)   -  nodeLocationTable(element(j),1));
            
            % %Get the coefficient for a_ij locally 
            % a(i,l) = - k_Array(e) * (dNj_dx*dNi_dx + dNj_dy*dNi_dy)/(4 * A);

            a(i,l) = - k_air * (dNj_dx*dNi_dx + dNj_dy*dNi_dy)/(4 * A);

            %get global matrix by placing (summing) the current coeff at find(element==i),find(element==l)
            globalMatrix(element(i),element(l)) = globalMatrix(element(i),element(l)) + a(i,l);

        end 
    end 

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specifying boundary conditions
% Here If it is a boundary, I know the node will be a certain constant
% temperature (constant temperature boundary condition Type 1). This means the row can be set to [0...1...0]  * T_BC = Q
% OR 
% Apply the third type boundary condition Use the outside temeprature value

%% Dirichlet Boundary conditions for every wall: Working
% for i = 1: length(globalMatrix)
%     if boundary_Array(i) ~= 0 ]
%         globalMatrix(i,:) = [zeros(length(1:i)-1,1)', 1 , zeros(length(globalMatrix)-i,1)']; % Set constant temperature boundary (Type1)
%     end
% end
% rhs = boundary_Array; % define the RHS of the equation (Boundary conditions, fluxes)



%% Setting Leaking boundary for other walls 
%  Addition of heat flow at the nodes due to q_out through the dx
%  or dy from the heat differential between the outside temperature (known) and the
%  inside temperature (unknown)
%  Result - New terms to the Global matrix 
%           Add heat flux to RHS of equation based on known T_out

for i = 1:length(globalMatrix)
    % K*L/t * T_i 
    % arrayHeatFlux = k_node(i)/thick * 
    % We will get a 3 column vector : 
    %       K is conductivity of surface, 
    %       t is its thicknesss
    %       L1 is the distance between the node and previous one
    %       L2 is the distance between the node and next one
    current_X = nodeLocationTable(i,1);
    current_Y = nodeLocationTable(i,2);
        % We are not at the first node or the last node 
    if (current_X < .001) && (current_Y > .001) && (current_Y < y_max - .001)
        % We are on the y axis either at x = 0  Vertical
        % direction in cross section of the room
        % Also not in the corners
        
        %Index of nodes that are along the y axis:
        currentArrayIndex = find(nodeLocationTable(:,1) <= .01);  
        % find the index in the current array of indexes that has the
        % value of interest 'nodeLocationTable(i,2)':
        indexOfMatch = find(abs(nodeLocationTable(currentArrayIndex,2) - nodeLocationTable(i,2)) <= .01 ); 

        % Find the index of the nodes next to i in terms of the
        % 'nodeLocationArray'
        nodeBefore  = currentArrayIndex(indexOfMatch-1);
        currentNode = i;
        nodeAfter   = currentArrayIndex(indexOfMatch+1);

        L1 =  abs(nodeLocationTable(i,2) - nodeLocationTable(nodeBefore,2)); 
        L2 =  abs(nodeLocationTable(i,2) - nodeLocationTable(nodeAfter,2)); 

        heatFlux_known = k_node(i)/ thickness_window .* [ L1/6, (L1 + L2)/3 , L2/6 ];

        globalMatrix(i,nodeBefore)  = globalMatrix(i,nodeBefore)  - heatFlux_known(1);
        globalMatrix(i,currentNode) = globalMatrix(i,currentNode) - heatFlux_known(2);
        globalMatrix(i,nodeAfter)   = globalMatrix(i,nodeAfter)   - heatFlux_known(3);

        wallHeat_asFlux_Array(i) = k_node(i)/ thickness_window  * (1/2)* (L1+L2) * boundaryTemperature;

    elseif (current_X > x_max - .00001) && (current_Y > .001) && (current_Y < y_max - .001)
        % We are on the y axis either at x = x_max  Vertical
        % direction in cross section of the room
        % Also not in the corners

        currentArrayIndex = find(nodeLocationTable(:,1) >= x_max-.001);  %Index of nodes that are along the x = xmax
        indexOfMatch = find(abs(nodeLocationTable(currentArrayIndex,2) - nodeLocationTable(i,2)) <= .01 );

        nodeBefore  = currentArrayIndex(indexOfMatch-1);
        currentNode = i;
        nodeAfter   = currentArrayIndex(indexOfMatch+1);

        L1 =  abs(nodeLocationTable(i,2) - nodeLocationTable(nodeBefore,2)); 
        L2 =  abs(nodeLocationTable(i,2) - nodeLocationTable(nodeAfter,2)); 

        heatFlux_known = k_node(i)/ thickness_window .* [ L1/6, (L1 + L2)/3 , L2/6 ];

        globalMatrix(i,nodeBefore)  = globalMatrix(i,nodeBefore)  - heatFlux_known(1);
        globalMatrix(i,currentNode) = globalMatrix(i,currentNode) - heatFlux_known(2);
        globalMatrix(i,nodeAfter)   = globalMatrix(i,nodeAfter)   - heatFlux_known(3);

        wallHeat_asFlux_Array(i) = k_node(i)/ thickness_window  * (1/2)* (L1+L2) * boundaryTemperature;

    elseif (current_Y < .001) && (current_X > .001) && (current_X < x_max - .001)
        % We are on the x axis at y = 0  horizontal
        % direction in cross section of the room
        % Also not in the corners
        
        currentArrayIndex = find(nodeLocationTable(:,2) <= .001);
        indexOfMatch = find(abs(nodeLocationTable(currentArrayIndex,1) - nodeLocationTable(i,1)) <= .01 );

        nodeBefore  = currentArrayIndex(indexOfMatch-1);
        currentNode = i;
        nodeAfter   = currentArrayIndex(indexOfMatch+1);

        L1 =  abs(nodeLocationTable(i,1) - nodeLocationTable(nodeBefore,1)); 
        L2 =  abs(nodeLocationTable(i,1) - nodeLocationTable(nodeAfter,1)); 

        heatFlux_known = k_node(i)/ thickness_window .* [ L1/6, (L1 + L2)/3 , L2/6 ];

        globalMatrix(i,nodeBefore)  = globalMatrix(i,nodeBefore)  - heatFlux_known(1);
        globalMatrix(i,currentNode) = globalMatrix(i,currentNode) - heatFlux_known(2);
        globalMatrix(i,nodeAfter)   = globalMatrix(i,nodeAfter)   - heatFlux_known(3);

        wallHeat_asFlux_Array(i) = k_node(i)/ thickness_window * (1/2)* (L1+L2) * boundaryTemperature;
   
    elseif (current_X < .001) && (current_Y < .001)
        %If bottom Left Corner:

        currentArrayIndex_X = find(nodeLocationTable(:,1) <= .001);             % array of indexes when x = 0 
        indexOfMatch_X      = find(abs(nodeLocationTable(currentArrayIndex_X,2) <= .01 )); % index of array above when y = 0 

        currentArrayIndex_Y   = find(nodeLocationTable(:,2) <= .001);% array of indexes when y = 0 
        indexOfMatch_Y      = find(abs(nodeLocationTable(currentArrayIndex_Y,1) <= .01 ));% index of array above when x = 0 

        nodeBefore  = currentArrayIndex_X(indexOfMatch_X+1);
        currentNode = i;
        nodeAfter   = currentArrayIndex_Y(indexOfMatch_Y+1);

        L1 =  abs(nodeLocationTable(i,1) - nodeLocationTable(currentArrayIndex_X(indexOfMatch_X+1),1)); 
        L2 =  abs(nodeLocationTable(i,2) -  nodeLocationTable(currentArrayIndex_Y(indexOfMatch_Y+1),2)); 

        heatFlux_known = k_node(i)/ thickness_window .* [ L1/6, (L1 + L2)/3 , L2/6 ];

        globalMatrix(i,nodeBefore)  = globalMatrix(i,nodeBefore)  - heatFlux_known(1);
        globalMatrix(i,currentNode) = globalMatrix(i,currentNode) - heatFlux_known(2);
        globalMatrix(i,nodeAfter)   = globalMatrix(i,nodeAfter)   - heatFlux_known(3);

        wallHeat_asFlux_Array(i) = k_node(i)/ thickness_window * (1/2)* (L1+L2) * boundaryTemperature;

     elseif (current_X > x_max - .001) && (current_Y > y_max - .001)
        %Corner 2 X large, y large - ZERO - Wall Boundary
        currentArrayIndex_X = find(nodeLocationTable(:,1) > x_max - .001 );
        indexOfMatch_X      = find(abs(nodeLocationTable(currentArrayIndex,2) > y_max - .001 ));

        currentArrayIndex_Y   = find(nodeLocationTable(:,2)  > y_max - .001);
        indexOfMatch_Y        = find(abs(nodeLocationTable(currentArrayIndex_Y,1) > x_max - .001 ));

        nodeBefore  = currentArrayIndex_X(indexOfMatch_X);
        currentNode = i;
        nodeAfter   = currentArrayIndex_Y(indexOfMatch_Y);

        L1 =  0; 
        L2 =  0; 

        heatFlux_known = k_node(i)/ thickness_window .* [ L1/6, (L1 + L2)/3 , L2/6 ];

        globalMatrix(i,nodeBefore)  = globalMatrix(i,nodeBefore)  - heatFlux_known(1);
        globalMatrix(i,currentNode) = globalMatrix(i,currentNode) - heatFlux_known(2);
        globalMatrix(i,nodeAfter)   = globalMatrix(i,nodeAfter)   - heatFlux_known(3);

        wallHeat_asFlux_Array(i) = k_node(i)/ thickness_window * (1/2)* (L1+L2) * boundaryTemperature;

     elseif (current_X < .001 ) && (current_Y > y_max - .001)
        % Corner X small and Y large - ZERO - Wall Boundary
        currentArrayIndex_X = find(nodeLocationTable(:,1) < .001 );
        indexOfMatch_X      = find(abs(nodeLocationTable(currentArrayIndex,2) > y_max - .001 ));

        currentArrayIndex_Y   = find(nodeLocationTable(:,2)  > y_max - .001);
        indexOfMatch_Y        = find(abs(nodeLocationTable(currentArrayIndex_Y,1)  <= .01));

        nodeBefore  = currentArrayIndex_X(indexOfMatch_X-1);
        currentNode = i;
        nodeAfter   = currentArrayIndex_Y(indexOfMatch_Y+1);

        L1 = 0; 
        L2 = 0; 

        heatFlux_known = k_node(i)/ thickness_window .* [ L1/6, (L1 + L2)/3 , L2/6 ];

        globalMatrix(i,nodeBefore)  = globalMatrix(i,nodeBefore)  - heatFlux_known(1);
        globalMatrix(i,currentNode) = globalMatrix(i,currentNode) - heatFlux_known(2);
        globalMatrix(i,nodeAfter)   = globalMatrix(i,nodeAfter)   - heatFlux_known(3);

        wallHeat_asFlux_Array(i) = k_node(i)/ thickness_window * (1/2)* (L1+L2) * boundaryTemperature;

      elseif (current_X > x_max - .001 ) && (current_Y < .001)
        % Corner X large and Y small

        currentArrayIndex_X = find(nodeLocationTable(:,1) > x_max - .001 );                % Get array of max x nodes
        indexOfMatch_X      = find(abs(nodeLocationTable(currentArrayIndex_X,2)  <= .01 )); % index when y is small

        currentArrayIndex_Y   = find(nodeLocationTable(:,2)  < .001);
        indexOfMatch_Y        = find(abs(nodeLocationTable(currentArrayIndex_Y,1) > x_max - .001 ));

        nodeBefore  = currentArrayIndex_X(indexOfMatch_X+1);
        currentNode = i;
        nodeAfter   = currentArrayIndex_Y(indexOfMatch_Y-1);

        L1 =   abs(nodeLocationTable(i,1) - nodeLocationTable(currentArrayIndex_X(indexOfMatch_X+1),1)); 
        L2 =   abs(nodeLocationTable(i,1) - nodeLocationTable(currentArrayIndex_Y(indexOfMatch_Y-1),2)); 

        heatFlux_known = k_node(i)/ thickness_window .* [ L1/6, (L1 + L2)/3 , L2/6 ];

        globalMatrix(i,nodeBefore)  = globalMatrix(i,nodeBefore)  - heatFlux_known(1);
        globalMatrix(i,currentNode) = globalMatrix(i,currentNode) - heatFlux_known(2);
        globalMatrix(i,nodeAfter)   = globalMatrix(i,nodeAfter)   - heatFlux_known(3);
        
        wallHeat_asFlux_Array(i) = k_node(i)/ thickness_window * (1/2) * (L1+L2) * boundaryTemperature; % 
    end
end

%% Derichlet BC for House facing wall + Leaky boundary for the outside facing walls
for i = 1:length(globalMatrix)
    %Setting the house facing wall conditions In global Matrix 
    if boundary_Array_House(i) ~= 0 
        globalMatrix(i,:) = [zeros(length(1:i)-1,1)', 1 , zeros(length(globalMatrix)-i,1)']; % Set constant temperature boundary (Type1)
    end
end
% rhs = boundary_Array_House + wallHeat_asFlux_Array; % define the RHS of the equation (Boundary conditions, fluxes)



%% Fluxes at certain nodes - Adding Heaters (manually for this)
heaterArray = zeros(length(globalMatrix),1);
close all
hold_Q_Loc = -2*ones(5,2);

%% 1 Heater (centered)

Q_heatSource  = 40;
scale   = 1;

denom  = 4;
num2   = 1;

q_loc = [1*wallLength2/2 num2*wallLength1/denom]; %Placing a heater at the center of the room
heater_node_index = dsearchn(nodeLocationTable,q_loc);
heaterArray(heater_node_index) = Q_heatSource;
hold_Q_Loc(5,:) = q_loc;

rhs = boundary_Array_House - wallHeat_asFlux_Array - heaterArray;

%% 3 Heaters
% 
% denom  = 10;
% num1   = 1;
% num2   = 9;
% Q_heatSource=10;
% middle_Q    = 1.2*Q_heatSource;
% scale       = 1;
% 
% q_loc = [num1*wallLength2/denom num1*wallLength1/denom]; %Placing a heater at the center of the room
% hold_Q_Loc(1,:) = q_loc;
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = scale*Q_heatSource;
% q_loc = [num2*wallLength2/denom num1*wallLength1/denom]; %Placing a heater at the center of the room
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = scale*Q_heatSource;
% hold_Q_Loc(2,:) = q_loc;
% 
% denom  = 2;
% num2   = 1;
% 
% q_loc = [num2*wallLength2/denom num2*wallLength1/denom]; %Placing a heater at the center of the room
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = middle_Q;
% hold_Q_Loc(5,:) = q_loc;
% rhs = boundary_Array_House - wallHeat_asFlux_Array - heaterArray;



%% 5 Heaters (4 corners and a middle)

% denom  = 10;
% num1   = 1;
% num2   = 9;
% Q_heatSource=5;
% middle_Q = 1.5*Q_heatSource;
% scale = 1;
% 
% hold_Q_Loc = zeros(5,2);
% 
% q_loc = [num1*wallLength2/denom num1*wallLength1/denom]; %
% hold_Q_Loc(1,:) = q_loc;
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = scale*middle_Q;
% 
% q_loc = [num2*wallLength2/denom num1*wallLength1/denom]; %
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = scale*middle_Q;
% hold_Q_Loc(2,:) = q_loc;
% 
% q_loc = [num1*wallLength2/denom num2*wallLength1/denom]; %
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = Q_heatSource;
% hold_Q_Loc(3,:) = q_loc;
% 
% q_loc = [num2*wallLength2/denom num2*wallLength1/denom]; %
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = Q_heatSource;
% hold_Q_Loc(4,:) = q_loc;
% 
% denom  = 2;
% num2   = 1;
% 
% q_loc = [num2*wallLength2/denom num2*wallLength1/denom]; %Placing a heater at the center of the room
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = middle_Q;
% hold_Q_Loc(5,:) = q_loc;
% rhs = boundary_Array_House - wallHeat_asFlux_Array - heaterArray;


%% 5 Heaters along walls

% denom  = 7;
% num1   = 1;
% num2   = 6;
% Q_heatSource=5;
% middle_Q = 1*Q_heatSource;
% scale = 1;
% 
% hold_Q_Loc = zeros(5,2);
% q_loc = [num1*wallLength2/denom num1*wallLength1/denom]; %Placing a heater at the center of the room
% hold_Q_Loc(1,:) = q_loc;
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = scale*Q_heatSource;
% q_loc = [num2*wallLength2/denom num1*wallLength1/denom]; %Placing a heater at the center of the room
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = scale*Q_heatSource;
% hold_Q_Loc(2,:) = q_loc;
% 
% q_loc = [num1*wallLength2/denom 4.5*wallLength1/denom]; %Placing a heater at the center of the room
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = Q_heatSource;
% hold_Q_Loc(3,:) = q_loc;
% q_loc = [num2*wallLength2/denom 4.5*wallLength1/denom]; %Placing a heater at the center of the room
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = Q_heatSource;
% hold_Q_Loc(4,:) = q_loc;
% 
% denom  = 2;
% num2   = 1;
% 
% q_loc = [num2*wallLength2/denom 1*wallLength1/4]; %Placing a heater at the center of the room
% heater_node_index = dsearchn(nodeLocationTable,q_loc);
% heaterArray(heater_node_index) = middle_Q;
% hold_Q_Loc(5,:) = q_loc
% rhs = boundary_Array_House - wallHeat_asFlux_Array - heaterArray;



%% Solving for temperature 
% invGLob = inv(globalMatrix);
temperature_C = inv(globalMatrix) * rhs; % Temperature solve for deg C 

%% % Display the results as a contour plot
close all
figure;
subplot(1,2,1);
trisurf(transformationTable,nodeLocationTable(:,1), nodeLocationTable(:,2), temperature_C);
title('Temperature Distribution');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Temperature (C)');
clim([-23 30])
% colorbar;
zlim([min(temperature_C) max(temperature_C)])

subplot(1,2,2);
% figure
trisurf(transformationTable,nodeLocationTable(:,1), nodeLocationTable(:,2), temperature_C);
title('Temperature Distribution With Heater Locations');
xlabel('X (m)');
ylabel('Y (m)');
colorbar;
clim([-23 30])
zlim([min(temperature_C) max(temperature_C)])
valHold  =max(temperature_C)
fondSize =40
text(hold_Q_Loc(1,1),hold_Q_Loc(1,2),valHold,". Q1",'Color','r','FontSize',fondSize,'FontWeight','bold')
text(hold_Q_Loc(2,1),hold_Q_Loc(2,2),valHold,". Q2",'Color','r','FontSize',fondSize,'FontWeight','bold')
text(hold_Q_Loc(3,1),hold_Q_Loc(3,2),valHold,". Q3",'Color','r','FontSize',fondSize,'FontWeight','bold')
text(hold_Q_Loc(4,1),hold_Q_Loc(4,2),valHold,". Q4",'Color','r','FontSize',fondSize,'FontWeight','bold')
text(hold_Q_Loc(5,1),hold_Q_Loc(5,2),valHold,". Q5",'Color','r','FontSize',fondSize,'FontWeight','bold')
fontsize(25,"points");
view(2)
hold on 
bufferSize = 10/39.37% Plotting a rectangle to represent the area we are trying to keep warm 20C
x_heated   = [bufferSize         4.8006-bufferSize    4.8006-bufferSize         bufferSize               bufferSize]
y_heated   = [ bufferSize         bufferSize          3.5814-bufferSize         3.5814-bufferSize        bufferSize]
plot3(x_heated,y_heated,valHold*ones(1,length(y_heated)), 'g', 'LineWidth', 5);


%% Average temperature of room
mean(temperature_C)
%% plotting heater locations:
% 
% figure;
% plot(x_rect, y_rect, 'b', 'LineWidth', 2);
% axis equal; % Make the aspect ratio equal
% title('Cross Section of Room');
% xlabel('X (m)');
% ylabel('Y (m)');
% grid on;
% hold on 
% bufferSize = 10/39.37 %10in from the wall /39.37 - > m
% x_heated = [bufferSize    4.8006-bufferSize    4.8006-bufferSize         bufferSize         bufferSize]
% y_heated = [ bufferSize         bufferSize    3.5814-bufferSize    3.5814-bufferSize        bufferSize]
% plot(x_heated,y_heated, 'g', 'LineWidth', 2);
% hold on 
% text(hold_Q_Loc(1,1),hold_Q_Loc(1,2),". Q1",'Color','r','FontSize',fondSize,'FontWeight','bold')
% text(hold_Q_Loc(2,1),hold_Q_Loc(2,2),". Q2",'Color','r','FontSize',fondSize,'FontWeight','bold')
% text(hold_Q_Loc(3,1),hold_Q_Loc(3,2),". Q3",'Color','r','FontSize',fondSize,'FontWeight','bold')
% text(hold_Q_Loc(4,1),hold_Q_Loc(4,2),". Q4",'Color','r','FontSize',fondSize,'FontWeight','bold')
% text(hold_Q_Loc(5,1),hold_Q_Loc(5,2),". Q5",'Color','r','FontSize',fondSize,'FontWeight','bold')
% fontsize(25,"points");

