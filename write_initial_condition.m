function f = write_initial_condition(hu,hv,Z,Zx,Zy,Zcos,u_inlet,v_inlet)

ndof   = length(hu(:,1));
nelems = length(hu(1,:));
nnodes = length(Z);

fid = fopen('DI_lava_ic.txt','w');
fprintf(fid,'%s\n','initial conditions');
fprintf(fid,'%-8.0f %8.0f %8.0f\n',[ndof nelems nnodes]);

for i = 1:nelems
    for j = 1:ndof
        fprintf(fid,'%-8.0f %16.10e %16.10e\n',...
                [i hu(j,i) hv(j,i)]);
        
    end
end

for i = 1:nnodes
    fprintf(fid,'%-8.0f %16.10e  %16.10e %16.10e %16.10e\n', ...
                [i Z(i) Zx(i) Zy(i) Zcos(i)]);
end

for i = 1:nnodes
    fprintf(fid,'%-8.0f %16.10e %16.10e\n',[i  u_inlet(i) v_inlet(i)]);
end

f = fclose('all');