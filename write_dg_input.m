function f = write_dg_input(A,B,C,Sc,PHI,PHInode,PHIedge,PSI,PSIedge,p)

ndof   = length(A(:,1,1));
napts  = length(A(1,:,1));
nlpts  = length(B(1,:,1));
nspts  = length(Sc(:,1));

fid = fopen('dg_input.txt','w');
fprintf(fid,'%s\n','dg input');
fprintf(fid,'%-8.0f %8.0f %8.0f %8.0f %8.0f\n',[ndof napts nlpts nspts p]);

for i = 1:2
    for j = 1:napts
        for k = 1:ndof
            fprintf(fid,'%16.10e\n', A(k,j,i));
        end
    end
end

for i = 1:3
    for j = 1:nlpts
        for k = 1:ndof
            fprintf(fid,'%16.10e\n', B(k,j,i));
        end
    end
end


for i = 1:napts
    for j = 1:ndof
        fprintf(fid,'%16.10e\n', C(j,i));
    end
end

for i = 1:nspts
    for j = 1:ndof
        fprintf(fid,'%16.10e\n', Sc(j,i));
    end
end

for i = 1:ndof
    for j = 1:napts
        fprintf(fid,'%16.10e\n', PHI(j,i));
    end
end

for i = 1:ndof
    for j = 1:3
        fprintf(fid,'%16.10e\n', PHInode(j,i));
    end
end

for i = 1:3
    for j = 1:ndof
        for k = 1:nlpts
            fprintf(fid,'%16.10e\n', PHIedge(k,j,i));
        end
    end
end

for i = 1:3
    for j = 1:napts 
        fprintf(fid,'%16.10e\n', PSI(j,i));
    end
end

for i = 1:3
    for j = 1:3
        for k = 1:nlpts
            fprintf(fid,'%16.10e\n', PSIedge(k,j,i));
        end
    end
end




f = fclose('all');