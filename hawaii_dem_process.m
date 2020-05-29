load dem_hawaii.txt
y_shift = max(YNODES)-max(dem_hawaii(:,2));
x_shift = 0;
dem_hawaii_xyz = zeros(size(dem_hawaii));
dem_hawaii_xyz(:,1) = dem_hawaii(:,1) - min(dem_hawaii(:,1))-x_shift;
dem_hawaii_xyz(:,2) = dem_hawaii(:,2) + y_shift;%+ 1.3466e+04;
dem_hawaii_xyz(:,3) = dem_hawaii(:,3); 
z_min = min(dem_hawaii_xyz(:,3));
z_shift = z_min - hint;
dem_hawaii_xyz(:,3) = dem_hawaii_xyz(:,3) - z_shift;
Xdem = reshape(dem_hawaii_xyz(:,1),39,38)';
Ydem = reshape(dem_hawaii_xyz(:,2),39,38)'; 
Zdem = reshape(dem_hawaii_xyz(:,3),39,38)'; 
ZX   = zeros(size(Xdem));
ZY   = zeros(size(Ydem));
[LY,LX] = size(Xdem);
dx = (Xdem(1,2)-Xdem(1,1));
dy = (Ydem(2,1)-Ydem(1,1));
twodeltax = 1/(2*dx);
twodeltay = 1/(2*dy);
% Inner boundary
i = (2:LX-1); j = (2:LY-1);
ZX(j,i) = (Zdem(j,i+1)-Zdem(j,i-1))*twodeltax;
ZY(j,i) = (Zdem(j+1,i)-Zdem(j-1,i))*twodeltay;

% Left boundary
i = 1; j = (2:LY-1);
ZX(j,i) = (-3*Zdem(j,i)+4*Zdem(j,i+1)-Zdem(j,i+2))*twodeltax;
ZY(j,i) = (Zdem(j+1,i)-Zdem(j-1,i))*twodeltay;

% Right boundary
i = LX; j = 2:LY-1;
ZX(j,i) = (Zdem(j,i-2)-4*Zdem(j,i-1)+3*Zdem(j,i))*twodeltax;
ZY(j,i) = (Zdem(j+1,i)-Zdem(j-1,i))*twodeltay;


% Top boundary
j = LY; i = (2:LX-1);
ZX(j,i) = (Zdem(j,i+1)-Zdem(j,i-1))*twodeltax;
ZY(j,i) = (Zdem(j-2,i)-4*Zdem(j-1,i)+3*Zdem(j,i))*twodeltay;


% Bottom boundary
j = 1; i = (2:LX-1);
ZX(j,i) = (Zdem(j,i+1)-Zdem(j,i-1))*twodeltax;
ZY(j,i) = (-3*Zdem(j,i)+4*Zdem(j+1,i)-Zdem(j+2,i))*twodeltay;


% Top left corner
i = 1;
j = LY;
ZX(j,i) = (-3*Zdem(j,i)+4*Zdem(j,i+1)-Zdem(j,i+2))*twodeltax;
ZY(j,i) = (Zdem(j-2,i)-4*Zdem(j-1,i)+3*Zdem(j,i))*twodeltay;

% Bottom left corner
i = 1;
j = 1;
ZX(j,i) = (-3*Zdem(j,i)+4*Zdem(i+1,j)-Zdem(i+2,j))*twodeltax;
ZY(j,i) = (-3*Zdem(j,i)+4*Zdem(j+1,i)-Zdem(j+2,i))*twodeltay;

% Bottom right corner
i = LX;
j = 1;
ZX(j,i) = (Zdem(j,i-2)-4*Zdem(j,i-1)+3*Zdem(j,i))*twodeltax;
ZY(j,i) = (-3*Zdem(j,i)+4*Zdem(j+1,i)-Zdem(j+2,i))*twodeltay;

% Top right corner
i = LX;
j = LY;
ZX(j,i) = (Zdem(j,i-2)-4*Zdem(j,i-1)+3*Zdem(j,i))*twodeltax;
ZY(j,i) = (Zdem(j-2,i)-4*Zdem(j-1,i)+3*Zdem(j,i))*twodeltay;


zX = interp2(Xdem,Ydem,ZX,XNODES,YNODES,'linear',-0.03);
zY = interp2(Xdem,Ydem,ZY,XNODES,YNODES,'linear',-0.03);
zD = interp2(Xdem,Ydem,Zdem,XNODES,YNODES);

fid = fopen('dem_hawaii_xyz.xyz','w');


for i = 1:length(dem_hawaii_xyz(:,1))

        fprintf(fid,'%16.10e %16.10e %16.10e\n',...
                [Xdem(:,1) Ydem(:,2) Zdem(:,3)]);
        
end

f = fclose('all');
% xdem = sort(dem_hawaii_xyz(:,1));
% ydem = sort(dem_hawaii_xyz(:,2));
% [Xdem,Ydem] = meshgrid(xdem,ydem);
% Zdem = interp2(dem_hawaii_xyz(:,1),dem_hawaii_xyz(:,2),dem_hawaii_xyz(:,3),...
%        Xdem,Ydem);
% mesh(Xdem,Ydem,Zdem)

    