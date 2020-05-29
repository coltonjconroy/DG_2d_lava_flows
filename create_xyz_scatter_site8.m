% interpolates Hawaii channel data to a background xyz grid 
load 'hawaii_site8_fine_bpts.mat'
offset = 111;
offset1 = 111;
avg_slope_1 = 2.8;
avg_slope_2 = 3.2;
nx = 250;
ny = 250;
xi = min(XNODES)-1; xf = max(XNODES)+1;
yi = 0 ;            yf = max(YNODES);
dx = (xf-xi)/nx;  dy = (yf-yi)/ny;
Xm = xi:dx:xf;    Ym = yi:dy:yf;
[x,y] = meshgrid(Xm,Ym);
% constraints for elevation 
% main channel cut-off
yM1 = main_pt1.Position(2);
xM1 = main_pt1.Position(1);
yM2 = main_pt2.Position(2);
xM2 = main_pt2.Position(1);
sM  = (yM2 - yM1)/(xM2 - xM1);
% channel split cut-off
yS1 = split_pt1.Position(2);
xS1 = split_pt1.Position(1);
yS2 = split_pt2.Position(2);
xS2 = split_pt2.Position(1);
sS  = (yS2 - yS1)/(xS2 - xS1);
% flow center-lines
% left branch center-line
y_lb1 = left_pt1.Position(2);
x_lb1 = left_pt1.Position(1);
y_lb2 = left_pt2.Position(2);
x_lb2 = left_pt2.Position(1);
s_lb  = (y_lb2 - y_lb1)/(x_lb2 - x_lb1);
% right branch center-line
y_rb1 = right_pt1.Position(2);
x_rb1 = right_pt1.Position(1);
y_rb2 = right_pt2.Position(2);
x_rb2 = right_pt2.Position(1);
s_rb  = (y_rb2 - y_rb1)/(x_rb2 - x_rb1);
% slope interpolation
fid = fopen('hawaii_slopeX.xyz','w');
fid1= fopen('hawaii_slopeY.xyz','w');
fid2= fopen('hawaii_dhdx.xyz','w');
fid3= fopen('hawaii_dhdy.xyz','w');
sZx = zeros(size(x)); cZ = ones(size(x));
sZy = zeros(size(y)); 
elev= zeros(size(x));
dist= zeros(size(x));
dhdx = zeros(size(x));
dhdy = zeros(size(y));
MAD  = 100*ones(size(x));
S    = zeros(size(x));
slope_b1 = zeros(length(branch_1(:,1)),1);
slope_b2 = zeros(length(branch_2(:,1)),1);
ref_slope1 = 0; ref_slope2 = 0;
for i = 1:length(branch_1(:,1))-1;
    if i > 1
        ref_slope1 = slope_b1(i-1);
        ref_slope2 = slope_b2(i-1);
    end
        DY1 = (branch_1(i+1,2)-branch_1(i,2));
        DX1 = (branch_1(i+1,1)-branch_1(i,1));
        DY2 = (branch_2(i+1,2)-branch_2(i,2));
        DX2 = (branch_2(i+1,1)-branch_2(i,1));
        slope_b1(i) = atand(DY1/DX1) ;
        slope_b2(i) = atand(DY2/DX2); 
end
k   = 0;
for i = 1:nx
    for j = 1:ny
        k = k + 1;
        check_1 = (x(j,i) - xM1)*sM + yM1;
        check_2 = (x(j,i) - xS1)*sS + yS1;
        if y(j,i) <= (check_1)
            sZx(j,i) = 0;
            sZy(j,i) = 0;
            elev(j,i) = max(branch_1(:,2));
            dhdx(j,i) = 0;
            dhdy(j,i) = 0;
            S(j,i)    = 0;
        else
            if y(j,i) > check_2 % left branch (branch 1)
                sZx(j,i) = sind(3.8)*cos(atan(s_lb));
                sZy(j,i) = sind(3.8)*sin(atan(s_lb));
                cZ(j,i)  = cosd(3.8);
                % centerline check
%                 check_l = (x(j,i)-x_lb1)*s_lb + y_lb1;
%                 if y(j,i) > check_l % pt is above centerline
%                     theta = atan(s_lb);
%                     dx1   = x(j,i) - x_lb1;
%                     d1    = dx1/cos(theta);
%                     dy1   = dx1*tan(theta);
%                     dy2   = (y(j,i)-y_lb1) - dy1;
%                     d2    = dy2*sin(theta);
%                     mad   = dy2*cos(theta);
%                     dc    = d1 + d2 +offset1; 
%                     if dc < 0
%                         dc = 0;
%                     elseif dc > max(branch_1(:,1))
%                         dc = max(branch_1(:,1));
%                     end
%                     dist(j,i) = dc;
%                     MAD(j,i)  = mad;
%                     elev(j,i) = interp1(branch_1(:,1),branch_1(:,2),dist(j,i),'linear');
%                     S(j,i)    = interp1(branch_1(:,1),slope_b1,dist(j,i),'linear');
%                     if i >1 && j>1 && i < nx && j < ny
%                         dhdx(j,i) = elev(j,i)-elev(j,i-1);
%                         dhdy(j,i) = elev(j,i)-elev(j-1,i);
%                         DdDx      = dist(j,i)-dist(j,i-1);
%                         DdDy      = dist(j,i)-dist(j-1,i);
%                         if abs(DdDx) < 1e-6
%                             theta_x = 0;
%                         else
%                             theta_x   = atand(dhdx(j,i)/DdDx);
%                         end
%                         if abs(DdDy) < 1e-6
%                             theta_y = 0;
%                         else
%                             theta_y   = atand(dhdy(j,i)/DdDy);
%                         end
%                          sZx(j,i) = sind(S(j,i))*cos(theta);
%                          sZy(j,i) = sind(S(j,i))*sin(theta);
%                          cZ(j,i) = cosd(S(j,i));
%                     end
%                 else % pt is below centerline
%                     %keyboard
%                     if y(j,i) < y_lb1
%                         theta = atan(s_lb);
%                         dx1   = x(j,i) - x_lb1;
%                         dy2   = y_lb1 - y(j,i);
%                         dy1   = dx1*tan(theta);
%                         d2    = (dy1 + dy2)*sin(theta);
%                         mad   = (dy1 + dy2)*cos(theta);
%                         dd    = sqrt((x(j,i)-x_lb1)^2+(dy1)^2);
%                         dc    = dd - d2 +offset1;
%                         if dc < 0
%                             dc = 0;
%                         elseif dc > max(branch_1(:,1))
%                             dc = max(branch_1(:,1));
%                         end
%                         dist(j,i) = dc ;
%                     else
%                         theta = atan(s_lb);
%                         dx1   = (x(j,i)-x_lb1);
%                         yt    = dx1*tan(theta);
%                         dt    = sqrt(dx1^2+yt^2);
%                         mad   = (yt-(y(j,i)-y_lb1))*cos(theta);
%                         dc    = dt - (yt-(y(j,i)-y_lb1))*sin(theta) +offset1;
%                         if dc < 0
%                             dc = 0;
%                         elseif dc > max(branch_1(:,1))
%                             dc = max(branch_1(:,1));
%                         end
%                         dist(j,i) = dc ;
%                         MAD(j,i)  = mad;
%                     end
%                     elev(j,i) = interp1(branch_1(:,1),branch_1(:,2),dist(j,i),'linear');
%                     S(j,i)    = interp1(branch_1(:,1),slope_b1,dist(j,i),'linear');
%                     if i >1 && j>1 && i < nx && j < ny
%                         dhdx(j,i) = elev(j,i)-elev(j,i-1);
%                         dhdy(j,i) = elev(j,i)-elev(j-1,i);
%                         DdDx      = dist(j,i)-dist(j,i-1);
%                         DdDy      = dist(j,i)-dist(j-1,i);
%                         if abs(DdDx) < 1e-6
%                             theta_x = 0;
%                         else
%                             theta_x   = atand(dhdx(j,i)/DdDx);
%                         end
%                         if abs(DdDy) < 1e-6
%                             theta_y = 0;
%                         else
%                             theta_y   = atand(dhdy(j,i)/DdDy);
%                         end
%                         sZx(j,i) = sind(S(j,i))*cos(theta);
%                         sZy(j,i) = sind(S(j,i))*sin(theta);
%                         cZ(j,i) = cosd(S(j,i));
%                     end
%                 end
            else % right branch (branch 2)
%                sZx(j,i) = sind(3.2)*cos(atan(s_rb));
%                sZy(j,i) = sind(3.2)*sin(atan(s_rb));
                elev(j,i) = max(branch_2(:,2));
                dhdx(j,i) = 0;
                dhdy(j,i) = 0;
                % centerline check
                check_r = (x(j,i)-x_rb1)*s_rb + y_rb1;
                if y(j,i) > check_r % pt is above centerline
                    theta = atan(s_rb);
                    dx1   = x(j,i) - x_rb1;
                    d1    = dx1/cos(theta);
                    dy1   = dx1*tan(theta);
                    dy2   = (y(j,i)-y_rb1) - dy1;
                    d2    = dy2*sin(theta);
                    mad   = dy2*cos(theta);
                    dc    = d1 + d2 +offset; 
                    if dc < 0
                        dc = 0;
                    elseif dc > max(branch_2(:,1))
                        dc = max(branch_2(:,1));
                    end
                    dist(j,i) = dc;
                    MAD(j,i)  = mad;
                    elev(j,i) = interp1(branch_2(:,1),branch_2(:,2),dist(j,i),'linear');
                    S(j,i)    = interp1(branch_2(:,1),slope_b2,dist(j,i),'linear');
                    if i >1 && j>1 && i < nx && j < ny
                        dhdx(j,i) = elev(j,i)-elev(j,i-1);
                        dhdy(j,i) = elev(j,i)-elev(j-1,i);
                        DdDx      = dist(j,i)-dist(j,i-1);
                        DdDy      = dist(j,i)-dist(j-1,i);
                        theta_x   = atand(dhdx(j,i)/DdDx);
                        theta_y   = atand(dhdy(j,i)/DdDy);
                        sZx(j,i) = sind(S(j,i))*cos(theta);
                        sZy(j,i) = sind(S(j,i))*sin(theta);
                        cZ(j,i) = cosd(S(j,i));
                    end
                else % pt is below centerline
                    if y(j,i) < y_rb1
                        theta = atan(s_rb);
                        dy2   = y_rb1 - y(j,i);
                        dx1   = (x(j,i)-x_rb1);
                        dy1   = dx1*tan(theta);
                        d2    = (dy1 + dy2)*sin(theta);
                        mad   = (dy1 + dy2)*cos(theta);
                        dd    = sqrt(dx1^2+dy1^2);
                        dc    = dd - d2 +offset;
                        if dc < 0
                            dc = 0;
                        elseif dc > max(branch_2(:,1))
                            dc = max(branch_2(:,1));
                        end
                        dist(j,i) = dc;
                        MAD(j,i)  = mad;
                    else
                        theta = atan(s_rb);
                        dx1   = (x(j,i)-x_rb1);
                        dy2   = (y(j,i)-y_rb1);
                        yt    = dx1*tan(theta);
                        mad   = (yt-dy2)*cos(theta);
                        dt    = sqrt(dx1^2+yt^2);
                        dc    = dt - (yt-dy2)*sin(theta)+offset;
                        if dc < 0
                            dc = 0;
                        elseif dc > max(branch_2(:,1))
                            dc = max(branch_2(:,1));
                        end
                        dist(j,i) = dc;
                        MAD(j,i)  = mad;
                    end
                    elev(j,i) = interp1(branch_2(:,1),branch_2(:,2),dist(j,i),'linear');
                    S(j,i)    = interp1(branch_2(:,1),slope_b2,dist(j,i),'linear');
                    if i >1 && j>1 && i < nx && j < ny
                        dhdx(j,i) = elev(j,i)-elev(j,i-1);
                        dhdy(j,i) = elev(j,i)-elev(j-1,i);
                        DdDx      = dist(j,i)-dist(j,i-1);
                        DdDy      = dist(j,i)-dist(j-1,i);
                        theta_x   = atand(dhdx(j,i)/DdDx);
                        theta_y   = atand(dhdy(j,i)/DdDy);
                        sZx(j,i) = sind(S(j,i))*cos(theta);
                        sZy(j,i) = sind(S(j,i))*sin(theta); 
                        cZ(j,i) = cosd(S(j,i));
                    end
                end
            end
        end
        fprintf(fid, '%16.10e %16.10e %16.10e\n', [x(j,i) y(j,i) sZx(j,i)]);
        fprintf(fid1,'%16.10e %16.10e %16.10e\n', [x(j,i) y(j,i) sZy(j,i)]);
        fprintf(fid2,'%16.10e %16.10e %16.10e\n', [x(j,i) y(j,i) dhdx(j,i)]);
        fprintf(fid3,'%16.10e %16.10e %16.10e\n', [x(j,i) y(j,i) dhdy(j,i)]);
    end
end
f = fclose('all');