fid = fopen('parameter_input.txt','w');
b_layer = 1; 
alph    = n; 
b0      = 2.73e4;      % variables no longer used in main code. Need to update. 
mu0     = 6.025e-7;
fprintf(fid,'%s\n','parameter input');
fprintf(fid,'%8.0f\n',wall_bc_type);
fprintf(fid,'%16.10e\n',rho);
fprintf(fid,'%16.10e\n',hint);
fprintf(fid,'%16.10e\n',b_layer);
fprintf(fid,'%16.10e\n',alph);
fprintf(fid,'%16.10e\n',A);
fprintf(fid,'%16.10e\n',B);
fprintf(fid,'%16.10e\n',C);
fprintf(fid,'%16.10e\n',tau_yield);
fprintf(fid,'%16.10e\n',T_int);
fprintf(fid,'%16.10e\n',T_wall);
fprintf(fid,'%16.10e\n',T_crust);
fprintf(fid,'%16.10e\n',T_air);
fprintf(fid,'%16.10e\n',eff);
fprintf(fid,'%16.10e\n',b0);
fprintf(fid,'%16.10e\n',mu0);
fprintf(fid,'%16.10e\n',kt);
fprintf(fid,'%8.0f\n',vid);
fprintf(fid,'%8.0f\n',vid_frame_rate);
fprintf(fid,'%8.0f\n',limiter);

f = fclose('all');
