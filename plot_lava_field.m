figure
velocity = (u_bar_n.^2+v_bar_n.^2).^(1/2);
trisurf(CONN,XNODES,YNODES,velocity); view(2)
if strcmpi(v_field,'on')
    hold on; quiver3(XNODES,YNODES,(max(velocity)+1)*ones(size(XNODES)),u_bar_n,v_bar_n,zeros(size(XNODES)),'w')
end
caxis([Vmin Vmax]);
colorbar
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Velocity (m/s) and direction','FontSize',14,'Interpreter','Latex')
figure; trisurf(CONN,XNODES,YNODES,H_n,'Linestyle','None'); view(2)
colorbar
caxis([Hmin Hmax]);
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Thickness (m)','FontSize',14,'Interpreter','Latex')
figure; trisurf(CONN,XNODES,YNODES,t_bar_n); view(2)
colorbar
caxis([Tmin Tmax]);
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Temperature (K)','FontSize',14,'Interpreter','Latex')
figure; trisurf(CONN,XNODES,YNODES,mu_n,'Linestyle','None'); view(2)
colorbar
caxis([Mmin Mmax]);
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Viscosity ($Pa \cdot s$)','FontSize',14,'Interpreter','Latex')
figure; trisurf(CONN,XNODES,YNODES,dudy_n); view(2)
colorbar
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Derivative $\partial u / \partial y$ ','FontSize',14,'Interpreter','Latex')
figure; trisurf(CONN,XNODES,YNODES,dvdx_n); view(2)
colorbar
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Derivative $\partial v / \partial y$ ','FontSize',14,'Interpreter','Latex')
grad_h = sqrt(dhdx_n.^2+dhdy_n.^2);
figure; trisurf(CONN,XNODES,YNODES,grad_h); view(2)
colorbar
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Derivative $\partial h / \partial x$ ','FontSize',14,'Interpreter','Latex')
hold on; quiver3(XNODES,YNODES,(max(grad_h)+1)*ones(size(XNODES)),dhdx_n,dhdy_n,zeros(size(XNODES)),'w')
% figure; trisurf(CONN,XNODES,YNODES,bl_n); view(2)
% colorbar
% xlabel('x-coordinate (m)','FontSize',13,'Interpreter','Latex')
% ylabel('y-coordinate (m)','FontSize',13,'Interpreter','Latex')
% title('Boundary layer thickness (m)','FontSize',14,'Interpreter','Latex')
figure; trisurf(CONN,XNODES,YNODES,w_n); view(2)
colorbar
ylabel('Longitude ($^{\circ}$N)','FontSize',13,'Interpreter','Latex')
xlabel('Latitude ($^{\circ}$W)','FontSize',13,'Interpreter','Latex')
title('Vertical velocity (m/s)','FontSize',14,'Interpreter','Latex')