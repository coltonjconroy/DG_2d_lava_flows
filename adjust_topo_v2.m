% Randomly adjust topography from original DEM
for i = 1:nnodes
    xc = XNODES(i);
    yc = YNODES(i);
    check_1 = (xc - xM1)*sM + yM1;
    check_2 = (xc - xS1)*sS + yS1;
    if yc <= (check_1)
        Zy(i) = Zy(i)*1.0;
        Zx(i) = Zx(i)*1.0;
    else
        if yc > check_2 % left branch (branch 1)
            Zy(i) = zy_cent(i);
            Zx(i) = zx_cent(i);
        else % right branch (branch 2)
            Zy(i) = Zy(i)*1.0;
            Zx(i) = Zx(i)*1.0;
        end
    end
    
end