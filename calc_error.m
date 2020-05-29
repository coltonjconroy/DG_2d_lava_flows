% Calculate errors
lava_vid_data
V_vid_mesh = V(XNODES,YNODES);
V_h        = sqrt(u_bar_n.^2 + v_bar_n.^2);
V_dif      = V_h - V_vid_mesh; 
L_inf      = max(V_dif)
iL = find(L_inf == V_dif);
x_loc = XNODES(iL)
rms_error  = sqrt(sum(V_dif.^2)/nnodes)

convert_to_latlon

figure; hold on;
trisurf(CONN,XNODES,YNODES,V_vid_mesh);
view(2)
caxis([Vmin Vmax]);
colorbar
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Video Speed (m/s)','FontSize',14,'Interpreter','Latex')

figure; hold on;
trisurf(CONN,XNODES,YNODES,V_dif,'LineStyle','None');
view(2)
caxis([-6.0 6.0]);
colorbar
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Speed difference $||V_h - V_{video}||$ (m/s)','FontSize',14,'Interpreter','Latex')

%title('$n = 1.50$','FontSize',13,'Interpreter','Latex')