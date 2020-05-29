ic_multx = ones(size(XNODES));
ic_multy = ones(size(YNODES));
for i = 1:nnodes
    xc = XNODES(i);
    yc = YNODES(i);
    check_1 = (xc - xM1)*sM + yM1;
    check_2 = (xc - xS1)*sS + yS1;
    if yc <= (check_1)
        % nothing
    else
        if yc > check_2 % left branch (branch 1)
            Zy(i) = Zy(i)*1.5;
        else % right branch (branch 2)
            Zx(i) = Zx(i)*0.20;
            Zy(i) = Zy(i)*0.20;
            ic_multy(i) = 0.5;
            ic_multx(i) = 3;
        end
    end
end