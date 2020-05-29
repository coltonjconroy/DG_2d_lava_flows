%-----------------------------------------------------------------
% Pre-process program for unstructured depth-integrated lava model
%
% written by: Colton J. Conroy
%                  @ APAM
%                  8.23.18
%-----------------------------------------------------------------
clear;
set_up = 'mesh';
mesh_type = 'new';
edge_normal_fig = 'on';
new_dem = 'on';
modify_topo = 'off';
new_dem_adjust = 'off';
hint = 3.1;
hvel = 2;
mux  = 2400;
muy  = 2400;
g    = 9.81;
rho  = 2085;
vel_int = 6.8; % avg inlet vel.
v_amp = 2/3;   % factor to amplify inlet velocity
save_mesh_data = 'on';
% inlet limit boundary pt: right: ()
%                           left: ()                        
if strcmpi(set_up,'pts')
    % load geographical coordinate data
    load Hawaii.mat
    % convert to xy for calcs
    R = 6378206.4;    % Constants for lat/lon conversion
    CPPLON = PTS.cpplon*pi/180;
    CPPLAT = PTS.cpplat*pi/180;
    % Convert lat/lon to meters
    PTS.Poly.x = R*(pi/180*PTS.Poly.x - CPPLON)*cos(CPPLAT);
    PTS.Poly.y = R*(pi/180*PTS.Poly.y);
    neg_value  = find(PTS.Poly.x < 0);
    neg_check  = isempty(neg_value);
    if neg_check == 0
        shift = min(PTS.Poly.x);
    else
        shift = 0;
    end
    PTS.Poly.x = PTS.Poly.x - shift; % shift so all x values are positive.
    save('Hawaii_site8.mat','PTS','Settings','xyzFun');
else
    % Open fort.14 mesh file (read only)
    if strcmpi(mesh_type,'old')
        fid = fopen('hawaii_30_3m_clean.14','r');
    else
        fid = fopen('hawaii_site8_velboundary_max8_min1.14','r');
    end
    % Read in mesh name
    AGRID = fgetl(fid);
    % Read in the number of elements and the number of nodes
    Nelems = fscanf(fid,'%i',1);
    Nnodes = fscanf(fid,'%i',1);
    % Read in node numbers and coordinates
    XNODES = zeros(Nnodes,1);
    YNODES = zeros(Nnodes,1);
    Z = zeros(Nnodes,1); 
    for i = 1 : Nnodes
        NODES(i,:) = (fscanf(fid,'%g',4))';
        XNODES(i) = NODES(i,2);
        YNODES(i) = NODES(i,3);
        Z(i)      = -NODES(i,4);
    end
    % Read in the element connectivity table
    CONN = zeros(Nelems,3);
    for j = 1 : Nelems
        ELEMENTS = (fscanf(fid,'%i',5))';
        CONN(j,:) = ELEMENTS(3:5);
    end
    fclose(fid);
    % generate mesh info (edge normals, element areas, etc.)
    rep  = TriRep(CONN,XNODES,YNODES);
    c_xy  = incenters(rep);
    c_bry = cartToBary(rep,[1:length(CONN(:,1))]',c_xy);
    [EDGES,ELEMS,NODES] = DG_meshData(rep);
    nelems = length(rep.Triangulation);
    nnodes = numel(NODES);
    nedges = numel(EDGES);
    if strcmpi(edge_normal_fig,'on')
        plot_mesh_normals
    end
    %--------------------------------
    % interpolate slope data to mesh
    %--------------------------------
    if strcmpi(mesh_type,'old')
        create_xyz_scatter_hawaii
        load dist_xy
    else
        %create_xyz_scatter_site8
        create_xyz_scatter_hawaii
        load dist_xy_site8
    end
    DhDx = interp2(x,y,dhdx,XNODES,YNODES);
    DhDy = interp2(x,y,dhdy,XNODES,YNODES);
    Zx   = interp2(x,y,sZx,XNODES,YNODES);
    Zy   = interp2(x,y,sZy,XNODES,YNODES);
    zx_cent = Zx;
    zy_cent = Zy;
    Zcos = interp2(x,y,cZ,XNODES,YNODES);
    Dc   = interp2(x,y,dist,XNODES,YNODES);
    Dmad = interp2(x,y,MAD,XNODES,YNODES);
    Dt   = interp2(X,Y,D,XNODES,YNODES);
    inlet_speed
    if strcmpi(mesh_type,'new')
        inlet_speed_from_video
    end
    nx_avg = sum(nx_inlet(:,2))/length(nx_inlet(:,2));
    ny_avg = sum(ny_inlet(:,2))/length(ny_inlet(:,2));
    %--------
    % new dem
    %--------
    if strcmpi(new_dem,'on')
        hawaii_dem_process
        clear Zx Zy
        Zx = zX;
        Zy = zY;
        izx = find(Zx > 0);
        izy = find(Zy > 0);
        Zx(izx) = -Zx(izx);
        Zy(izy) = -Zy(izy);
    end
    if strcmpi(modify_topo,'on')
        adjust_topo
    end
    if strcmpi(new_dem_adjust,'on')
        hawaii_dem_process
        zx_cent = Zx;
        zy_cent = Zy; 
        clear Zx Zy
        Zx = zX;
        Zy = zY;
        %izx = find(Zx > 0);
        %izy = find(Zy > 0);
        %Zx(izx) = -Zx(izx);
        %Zy(izy) = -Zy(izy);
        adjust_topo_v2
    end
     
    %--------------------
    % polynomial approx.
    %--------------------
    p = 0;
    ndof = (p+1)*(p+2)/2;
    %------------------------------
    % set-up dg matricies and basis
    %------------------------------
    dg_basis_and_matricies
    %----------------------------------
    % initial condition & L2-projection
    %----------------------------------
 %   if strcmpi(new_dem,'off')
 %       u_node = 0.30+(rho*g*Zx)./(3*mux).*hint.^2;
 %       v_node = 5.00+(rho*g*Zy)./(3*muy).*hint.^2;
 %   else
        u_node = 0.30-(rho*g*Zx)./(3*mux).*hvel.^2;
        v_node = 5.00-(rho*g*Zy)./(3*muy).*hvel.^2;
        if strcmpi(modify_topo,'on')
            u_node = u_node.*ic_multx;
            v_node = v_node.*ic_multy;
        end
 %   end
    VEL    = (u_node.^2 + v_node.^2).^(1/2);
    hu = zeros(ndof,nelems); hv = zeros(ndof,nelems);
    for j = 1:nelems
        
        U  = [u_node(ELEMS(j).nodes(1)); ...
              u_node(ELEMS(j).nodes(2)); ...
              u_node(ELEMS(j).nodes(3))];
        
        V  = [v_node(ELEMS(j).nodes(1)); ...
              v_node(ELEMS(j).nodes(2)); ...
              v_node(ELEMS(j).nodes(3))];
        
        u_l2 = Psi * U;
        v_l2 = Psi * V;
                
        hu(:,j) = C * u_l2 * hint;  % momentum hu
        hv(:,j) = C * v_l2 * hint;
    end
    clear Z; Z = hint*ones(size(XNODES));
    %------------------
    % write input files
    %------------------
    write_dg_input_v2(A,B,C,Sc,P,Ps,Pn,PHIedge,Psi,PSIedge,p)
    write_mesh_input(ELEMnodes,ELEMarea,ELEMxy,EDGEnormals,...
        EDGEelems,EDGEtype,EDGElengths,ELEMedges,EDGEnodes,XNODES,YNODES,...
        inlet_edges,outlet_edges,c_bry,SRCEpts,EDGEpts)
    write_initial_condition(hu,hv,Z,Zx,Zy,Zcos,u_inlet_nodes,v_inlet_nodes)
    if strcmpi(save_mesh_data,'on')
        if strcmpi(mesh_type,'old')
            save('Hawaii_run_mesh_info.mat','CONN','XNODES','YNODES','EDGES','ELEMS','ELEMENTS','NODES','hint','D','Dc');
        else
            save('Hawaii_run_mesh_info_site8.mat','CONN','XNODES','YNODES','EDGES','ELEMS','ELEMENTS','NODES','hint','Dc');
        end
    end
end