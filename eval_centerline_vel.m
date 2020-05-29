%load 'Hawaii_run_mesh_info.mat'
yS1 = split_pt1.Position(2);
xS1 = split_pt1.Position(1);
yS2 = split_pt2a.Position(2);
xS2 = split_pt2b.Position(1);
sS  = (yS2 - yS1)/(xS2 - xS1);
n   = 0; m = 0; 
for i = 1:nnodes
    if Dc(i) < 2.0
        y_check = (XNODES(i) - xS1)*sS + yS1;
        if YNODES(i) > y_check % upper branch            
            local_elem = [NODES(i).elems]';
            nelem_l    = length(local_elem);
            total_area = 0;
            for j = 1:nelem_l
                k = local_elem(j);
                total_area = total_area + [ELEMS(k).area];
            end
            u_local = 0; v_local = 0; h_local = 0; t_local = 0;
            for j = 1:nelem_l
                k = local_elem(j);
                H_l     = hint + Eta(1,k);
                u_bar   = hu_bar(1,k)/H_l;
                v_bar   = hv_bar(1,k)/H_l;
                t_bar   = ht_bar(1,k)/H_l;
                u_local = u_local + u_bar*[ELEMS(k).area]/total_area;
                v_local = v_local + v_bar*[ELEMS(k).area]/total_area;
                h_local = h_local + H_l*[ELEMS(k).area]/total_area;
                t_local = t_local + t_bar*[ELEMS(k).area]/total_area;
            end
            vel = sqrt(u_local^2 + v_local^2);
            n   = n + 1;
            left_branch(n,1) = vel;
            left_branch(n,2) = h_local;
            left_branch(n,3) = t_local;
            left_branch(n,4) = D(i); 
        elseif YNODES(i) < y_check % lower branch
            local_elem = [NODES(i).elems]';
            nelem_l    = length(local_elem);
            total_area = 0;
            for j = 1:nelem_l
                k = local_elem(j);
                total_area = total_area + [ELEMS(k).area];
            end
            u_local = 0; v_local = 0; h_local = 0; t_local = 0;
            for j = 1:nelem_l
                k = local_elem(j);
                H_l     = hint + Eta(1,k);
                u_bar   = hu_bar(1,k)/H_l;
                v_bar   = hv_bar(1,k)/H_l;
                t_bar   = ht_bar(1,k)/H_l;
                u_local = u_local + u_bar*[ELEMS(k).area]/total_area;
                v_local = v_local + v_bar*[ELEMS(k).area]/total_area;
                h_local = h_local + H_l*[ELEMS(k).area]/total_area;
                t_local = t_local + t_bar*[ELEMS(k).area]/total_area;
            end
            vel = sqrt(u_local^2 + v_local^2);
            m   = m + 1;
            right_branch(m,1) = vel;
            right_branch(m,2) = h_local;
            right_branch(m,3) = t_local;
            right_branch(m,4) = D(i);
        end
    end
end
% sort by distance
[Dl,li] = sort(left_branch(:,4));
[Dr,ri] = sort(right_branch(:,4));
vel_l = left_branch(li,1);
vel_r = right_branch(ri,1);

            
        