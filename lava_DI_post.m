%--------------------------------------------------------------------------
% program: lava_DI_post.m
%
% written by: Colton J. Conroy
%
% Reads in degrees of freedom from FORTRAN90 program and interpolates
% output to mesh nodes for visualization.
%
%--------------------------------------------------------------------------
clear;
figs = 'on';  % turns figures on or off
errors = 'on'; % calc errors?
vid  = 'off'; % make a video?
v_field = 'on'; % display velocity field?
vid_type = 'velocity'; % variable to display in video
vid_frame_rate = 10.0d0; % video frame rate
vid_file_name = 'Hawaii_velocity_avgslope'; % video name
load Hawaii_run_mesh_info_site8
nframes = 30; % number of frames for video
g = 9.81;
dt = 0.05; % time step used
% plot colorbar range
Tmin = 1400; Hmin = 1.0;   Mmin = 25;   Vmin = 0;
Tmax = 1425; Hmax = 20.0;  Mmax = 900;  Vmax = 9.0;
% read global output file
read_dof
% loop over nodes and interpolate values
u_bar_n = zeros(size(XNODES)); v_bar_n = zeros(size(XNODES));
H_n     = zeros(size(XNODES)); t_bar_n = zeros(size(XNODES));
mu_n    = zeros(size(XNODES)); dudy_n  = zeros(size(XNODES));
dvdx_n  = zeros(size(XNODES)); dhdx_n  = zeros(size(XNODES));
dhdy_n  = zeros(size(YNODES)); bl_n    = zeros(size(XNODES));
w_n     = zeros(size(XNODES)); dwdx_n  = zeros(size(XNODES));
dwdy_n  = zeros(size(XNODES));
for i = 1:nnodes
    local_elem = [NODES(i).elems]';
    nelem_l    = length(local_elem);
    total_area = 0;
    for j = 1:nelem_l
        k = local_elem(j);
        total_area = total_area + [ELEMS(k).area];
    end
    u_local = 0; v_local  = 0; h_local = 0; t_local = 0; mu_local = 0;
    du_local= 0; dv_local = 0; dhx_local = 0; dhy_local = 0; b_local = 0;
    w_local = 0; wx_local = 0; wy_local = 0; 
    for j = 1:nelem_l
        k = local_elem(j);
        H_l      = Eta(1,k);
        u_bar    = hu_bar(1,k)/H_l;
        v_bar    = hv_bar(1,k)/H_l;
        t_bar    = ht_bar(1,k)/H_l;
        u_local  = u_local  + u_bar*[ELEMS(k).area]/total_area;
        v_local  = v_local  + v_bar*[ELEMS(k).area]/total_area;
        h_local  = h_local  + H_l*[ELEMS(k).area]/total_area;
        t_local  = t_local  + t_bar*[ELEMS(k).area]/total_area;
        mu_local = mu_local + mu(1,k)*[ELEMS(k).area]/total_area;
        du_local = du_local + dudy(1,k)*[ELEMS(k).area]/total_area;
        dv_local = dv_local + dudy(1,k)*[ELEMS(k).area]/total_area;
        dhx_local= dhx_local+ dhdx(1,k)*[ELEMS(k).area]/total_area;
        dhy_local= dhy_local+ dhdy(1,k)*[ELEMS(k).area]/total_area;
        b_local  = b_local  + b_thick(1,k)*[ELEMS(k).area]/total_area;
        w_local  = w_local  + w(1,k)*[ELEMS(k).area]/total_area;
        wx_local = wx_local + dwdx(1,k)*[ELEMS(k).area]/total_area;
        wy_local = wy_local + dwdy(1,k)*[ELEMS(k).area]/total_area;
    end
    u_bar_n(i) = u_local;  v_bar_n(i) = v_local;
    H_n(i)     = h_local;  t_bar_n(i) = t_local;
    mu_n(i)    = mu_local;  dudy_n(i) = du_local;
    dvdx_n(i)  = dv_local;  dhdx_n(i) = dhx_local;
    dhdy_n(i)  = dhy_local; bl_n(i)   = b_local; 
    w_n(i)     = w_local; dwdx_n(i)   = wx_local;
    dwdy_n(i)  = wy_local; 
end

if strcmpi(errors,'on')
    calc_error
end

if strcmpi(figs,'on')
    plot_lava_field
end

if strcmpi(vid,'on')
    writerObj = VideoWriter(vid_file_name);
    writerObj.Quality = 100;
    open(writerObj);
    video
    close(writerObj);
end

