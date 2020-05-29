fid = fopen('output.txt','r');
nelems = fscanf(fid,'%i',1);
ndof   = fscanf(fid,'%i',1);
nnodes = fscanf(fid,'%i',1);

Eta = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        Eta(j,i) = fscanf(fid,'%f',1);
    end
end

hu_bar = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        hu_bar(j,i) = fscanf(fid,'%f',1);
    end
end

hv_bar = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        hv_bar(j,i) = fscanf(fid,'%f',1);
    end
end

ht_bar = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        ht_bar(j,i) = fscanf(fid,'%f',1);
    end
end

mu = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        mu(j,i) = fscanf(fid,'%f',1);
    end
end

dvdx = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        dvdx(j,i) = fscanf(fid,'%f',1);
    end
end

dudy = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        dudy(j,i) = fscanf(fid,'%f',1);
    end
end

dhdx = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        dhdx(j,i) = fscanf(fid,'%f',1);
    end
end

dhdy = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        dhdy(j,i) = fscanf(fid,'%f',1);
    end
end

b_thick = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        b_thick(j,i) = fscanf(fid,'%f',1);
    end
end

w = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        w(j,i) = fscanf(fid,'%f',1);
    end
end

dwdx = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        dwdx(j,i) = fscanf(fid,'%f',1);
    end
end

dwdy = zeros(ndof,nelems);
for i = 1:nelems
    for j = 1:ndof
        dwdy(j,i) = fscanf(fid,'%f',1);
    end
end

fclose(fid);