% inlet velocities
load inletspeeds.mat
x1_inlet = inlet_pt1.Position(1); x2_inlet = inlet_pt2.Position(1);
y1_inlet = inlet_pt1.Position(2); y2_inlet = inlet_pt2.Position(2);
x1_l_out = L_outlet_pt1.Position(1); x2_l_out = L_outlet_pt2.Position(1);
y1_l_out = L_outlet_pt1.Position(2); y2_l_out = L_outlet_pt2.Position(2);
x1_r_out = R_outlet_pt1.Position(1); x2_r_out = R_outlet_pt2.Position(1);
y1_r_out = R_outlet_pt1.Position(2); y2_r_out = R_outlet_pt2.Position(2);
l_anchor = D7(1,1);
D7(:,1)  = D7(:,1) - l_anchor;
l_inlet  = sqrt((x2_inlet-x1_inlet)^2 + (y2_inlet-y1_inlet)^2);
dxi = (l_inlet-x1_inlet)/(length(D7(:,1))-1);
di = x1_inlet:dxi:l_inlet;
slope_inlet = (y2_inlet - y1_inlet)/(x2_inlet - x1_inlet); % note this is only approx
theta_inlet = atan(slope_inlet); 
x_inlet = x1_inlet + di.*cos(theta_inlet);
y_inlet = y1_inlet + di.*sin(theta_inlet);
u_inlet_nodes = zeros(size(XNODES));
v_inlet_nodes = zeros(size(XNODES));
outlet_edges = zeros(nedges,1);
inlet_edges  = zeros(nedges,1);
k = 0;
for j = 1:nedges
    
    if EDGES(j).type == 1 % boundary edge
        
        nx = EDGES(j).normal(1);
        ny = EDGES(j).normal(2);
        
        nodes = EDGES(j).nodes;
        x_node1 = XNODES(nodes(1));
        y_node1 = YNODES(nodes(1));
        x_node2 = XNODES(nodes(2));
        y_node2 = YNODES(nodes(2));
        
        x_centroid = (x_node1 + x_node2)/2;
        y_centroid = (y_node1 + y_node2)/2;
        
        if x_centroid >= x1_inlet && x_centroid <= x2_inlet
            if y_centroid <= y1_inlet &&  y_centroid >= y2_inlet
                if nx < 0 && ny < 0
                    if x_node1 == x1_inlet
                        inlet_vel1 = D7(1,2);
                        inlet_vel2 = interp1(x_inlet,D7(:,2),x_node2);
                    elseif x_node2 == x1_inlet
                        inlet_vel2 = D7(1,2);
                        inlet_vel1 = interp1(x_inlet,D7(:,2),x_node1);
                    elseif x_node2 == x2_inlet
                        inlet_vel2 = D7(end,2);
                        inlet_vel1 = interp1(x_inlet,D7(:,2),x_node1);
                    elseif x_node2 == x2_inlet
                        inlet_vel2 = D7(end,2);
                        inlet_vel1 = interp1(x_inlet,D7(:,2),x_node1);
                    else
                        inlet_vel1 = interp1(x_inlet,D7(:,2),x_node1);
                        inlet_vel2 = interp1(x_inlet,D7(:,2),x_node2);
                    end
                    c1 = isnan(inlet_vel1);
                    c2 = isnan(inlet_vel2);
                    if c1 == 1 || c2 == 1
                        keyboard
                    end
                    u_inlet_nodes(nodes(1)) = inlet_vel1*nx;
                    u_inlet_nodes(nodes(2)) = inlet_vel2*nx;
                    v_inlet_nodes(nodes(1)) = inlet_vel1*ny;
                    v_inlet_nodes(nodes(2)) = inlet_vel2*ny;
                    inlet_edges(j) = 1;
                    scatter(x_centroid,y_centroid)
                    k = k + 1;
                    nx_inlet(k,1) = j;
                    nx_inlet(k,2) = nx; 
                    ny_inlet(k,1) = j;
                    ny_inlet(k,2) = ny;
                end
            end
        end
        if x_centroid >= x1_l_out && x_centroid <= x2_l_out
            if y_centroid <= y1_l_out &&  y_centroid >= y2_l_out
                if nx > 0 && ny > 0
                    outlet_edges(j) = 1;
                    scatter(x_centroid,y_centroid)
                end
            end
        end
        if x_centroid >= x1_r_out && x_centroid <= x2_r_out
            if y_centroid <= y1_r_out &&  y_centroid >= y2_r_out
                if nx > 0 && ny > 0
                    outlet_edges(j) = 1;
                    scatter(x_centroid,y_centroid)
                end
            end
        end
    end
    
end