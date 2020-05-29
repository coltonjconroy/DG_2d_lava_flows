
function [EDGES,ELEMS,NODES] = DG_meshData(rep)
%------------------------------------------------------------------------------------------------------
%   function [EDGES,ELEMS,NODES] = DG_meshData(rep)
%
%   DG_mesh_data stores mesh data for rep, the triangulation created in solver_Donelan.m
%
%   EDGES is a structure containing:
%   EDGES(i).elems - returns the two neighboring elements of global edge i
%   EDGES(i).nodes - returns the two global nodes of global edge i
%   EDGES(i).normal - returns the outward unit normal for global edge i
%   EDGES(i).length - returns the length of global edge i
%   EDGES(i).type -  0 for internal edge or 1 for boundary edge
%
%   ELEMS is a structure containing:
%   ELEMS(i).area - area of simplex i
%   ELEMS(i).nodes - global node numbers for local nodes 1,2,3 of element i
%   ELEMS(i).edges - global edge numbers for local edges 1,2,3 of element i
%   ELEMS(i).nbors - neighboring simplices for local edges 1,2,3 of element i
%  
%   NODES is a structure containing:
%   NODES(i).x - x-coordinate for global node i
%   NODES(i).y - y-coordinate for global node i
%   NODES(i).elems - returns all the element connected to global node i
%------------------------------------------------------------------------------------------------------
%   Written by Angela Nappi, October 2012

Nelems = length(rep.Triangulation);
Nnodes = length(rep.X);
Nedges = length(edges(rep));

%---------------------%
% Store node data %
%---------------------%

for i = 1 : Nnodes
    NODES(i).x = rep.X(i,1);
    NODES(i).y = rep.X(i,2);
    NODES(i).elems = cell2mat(vertexAttachments(rep,i));
end

%---------------------%
% Store edge data %
%---------------------%

E = edges(rep); % edge between 2 vertex coords
bList = freeBoundary(rep); % edges on the boundary
j = find(bList(:,1) > bList(:,2)); % reorganize
bList(j,:) = bList(j,[2,1]);
z = [0 0 1];

for i = 1 : Nedges
    EDGES(i).nodes = E(i,:);
    search = EDGES(i).nodes;
    search_bounds = any(ismember(bList,search,'rows'));
    if search_bounds == 1
        EDGES(i).type = 1;
    else
        EDGES(i).type = 0;
    end
    distx = NODES(EDGES(i).nodes(1)).x - NODES(EDGES(i).nodes(2)).x; 
    disty = NODES(EDGES(i).nodes(1)).y - NODES(EDGES(i).nodes(2)).y;
    EDGES(i).length = sqrt(distx^2 + disty^2);
    EDGES(i).elems = cell2mat(edgeAttachments(rep,E(i,:)));
end

%------------------------%
% Store element data %
%------------------------%

for i = 1 : Nelems
    ELEMS(i).nodes = rep.Triangulation(i,:);        
    ELEMS(i).x21 = NODES(ELEMS(i).nodes(2)).x - NODES(ELEMS(i).nodes(1)).x;
    ELEMS(i).x13 = NODES(ELEMS(i).nodes(1)).x - NODES(ELEMS(i).nodes(3)).x;
    ELEMS(i).y31 = NODES(ELEMS(i).nodes(3)).y - NODES(ELEMS(i).nodes(1)).y;
    ELEMS(i).y12 = NODES(ELEMS(i).nodes(1)).y - NODES(ELEMS(i).nodes(2)).y;
    ELEMS(i).nbors = neighbors(rep,i);
    if any(isnan(ELEMS(i).nbors)) == 1
        ELEMS(i).nbors(isnan(ELEMS(i).nbors) == 1) = 0; % Replace NaN with 0
    end
    L = rep(i,:);
    local_edge1 = [L(2),L(3)];
    local_edge1 = sort(local_edge1);
    ELEMS(i).edges(1) = find(E(:,1) == local_edge1(1) & E(:,2) == local_edge1(2));
    L1 = EDGES(ELEMS(i).edges(1)).length;
    local_edge2 = [L(1),L(3)];
    local_edge2 = sort(local_edge2);
    ELEMS(i).edges(2) = find(E(:,1) == local_edge2(1) & E(:,2) == local_edge2(2));
    L2 = EDGES(ELEMS(i).edges(2)).length;
    local_edge3 = [L(1),L(2)];
    local_edge3 = sort(local_edge3);
    ELEMS(i).edges(3) = find(E(:,1) == local_edge3(1) & E(:,2) == local_edge3(2));
    L3 = EDGES(ELEMS(i).edges(3)).length;
    R = 0.5*(L1+L2+L3);
    ELEMS(i).area = sqrt(R*(R-L1)*(R-L2)*(R-L3)); % Heron's formula
end

%------------------------%
% Find edge normals %
%------------------------%

for i = 1 : Nedges
    if EDGES(i).type == 0
        nx = NODES(EDGES(i).nodes(1)).x - NODES(EDGES(i).nodes(2)).x; 
        ny = NODES(EDGES(i).nodes(1)).y - NODES(EDGES(i).nodes(2)).y;
    elseif EDGES(i).type == 1
        nodelist = rep.Triangulation(EDGES(i).elems,:); 
        j = find(ELEMS(EDGES(i).elems).edges == i);
        if j == 1
            nx = NODES(nodelist(3)).x - NODES(nodelist(2)).x; 
            ny = NODES(nodelist(3)).y - NODES(nodelist(2)).y;
        elseif j == 2
            nx = NODES(nodelist(1)).x - NODES(nodelist(3)).x; 
            ny = NODES(nodelist(1)).y - NODES(nodelist(3)).y;
        elseif j == 3
            nx = NODES(nodelist(2)).x - NODES(nodelist(1)).x; 
            ny = NODES(nodelist(2)).y - NODES(nodelist(1)).y;
        end
    end
    nv = [nx ny 0]/EDGES(i).length; % normalize vector
    normal = cross(nv,z);
    EDGES(i).normal = normal(1,1:2);
end
  
end % function


