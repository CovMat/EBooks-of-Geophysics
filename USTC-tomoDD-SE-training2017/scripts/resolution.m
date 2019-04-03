%this script is used to calculate the resolvablity of the Lateral velocity
%resolution.

% the surrounding grid nodes are also used to calculate the Restoration degree of one node,
% thus you need to define the maximal distance from the surrounding nodes to the calculated node.
radius=0.5; 


%%%%%%%%%%% read velocity model %%%%%%%%%%%
vp=load('Vp_model.dat');
vs=load('Vs_model.dat');

%%%%%%%%%%%% read 1D MOD %%%%%%%%%%%%%%
fid=fopen('../MOD','r');
gridnumber=fscanf(fid,'%f',[1,4]);
nx=gridnumber(2);ny=gridnumber(3);nz=gridnumber(4);
X=fscanf(fid,'%f',[1,nx]);
Y=fscanf(fid,'%f',[1,ny]);
Z=fscanf(fid,'%f',[1,nz]);

for i=1:nx-1
    dx(i)=X(i+1)-X(i);
end
for i=1:ny-1
    dy(i)=Y(i+1)-Y(i);
end
for i=1:nz-1
    dz(i)=Z(i+1)-Z(i);
end

for k=1:nz
    for j=1:ny
        for i=1:nx
            vp_1d((k-1)*ny+j,i)=fscanf(fid,'%f',1);
        end
    end
end
for k=1:nz
    for j=1:ny
        for i=1:nx
            vpvs_1d((k-1)*ny+j,i)=fscanf(fid,'%f',1);
        end
    end
end
vs_1d=vp_1d./vpvs_1d;
fclose(fid);

%%%%%%%%%% read checkerbd MOD %%%%%%%%%%%%
fid=fopen('../Syn/MOD','r');
gridnumber=fscanf(fid,'%f',[1,4]);
nx=gridnumber(2);ny=gridnumber(3);nz=gridnumber(4);
X=fscanf(fid,'%f',[1,nx]);
Y=fscanf(fid,'%f',[1,ny]);
Z=fscanf(fid,'%f',[1,nz]);
for k=1:nz
    for j=1:ny
        for i=1:nx
            vp_checkbd((k-1)*ny+j,i)=fscanf(fid,'%f',1);
        end
    end
end
for k=1:nz
    for j=1:ny
        for i=1:nx
            vpvs_checkbd((k-1)*ny+j,i)=fscanf(fid,'%f',1);
        end
    end
end
vs_checkbd=vp_checkbd./vpvs_checkbd;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%% resolution P %%%%%%%%%%%%%%%%%%%%%%
res=ones(ny,nx,nz);

fid=fopen('res_P.dat','w');

for k=1:nz
    for i=1:nx
        for j=1:ny
            VEL(j,i,k)=vp((k-1)*ny+j,i);
            VEL_1d(j,i,k)=vp_1d((k-1)*ny+j,i);
            VEL_chk(j,i,k)=vp_checkbd((k-1)*ny+j,i);
        end
    end
end

for k=2:nz-1
    %orig=VEL(1,1,k);
    %disp(orig);
    for j=2:ny-1
        for i=2:nx-1
            mol=0.0;
            den=0.0;
           % distx=0;disty=0;
            for jj=2:ny-1
                for ii=2:nx-1
                    distx=0;disty=0;
                    if(i<ii)
                        ilow=i;iup=ii;
                    else
                        ilow=ii;iup=i;
                    end
                    if(j<jj)
                        jlow=j;jup=jj;
                    else
                        jlow=jj;jup=j;
                    end
                    for itmp=ilow:iup-1
                        distx=distx+dx(itmp);
                    end
                    for itmp=jlow:jup-1
                        disty=disty+dy(itmp);
                    end
                    dist=sqrt(distx^2+disty^2);
                    if(dist<=radius)
                        mol=mol+((VEL_chk(jj,ii,k)-VEL_1d(jj,ii,k))+(VEL(jj,ii,k)-VEL_1d(jj,ii,k)))^2;
                        den=den+(VEL_chk(jj,ii,k)-VEL_1d(jj,ii,k))^2+(VEL(jj,ii,k)-VEL_1d(jj,ii,k))^2;
                        %disp(mol);
                        %disp(den);
                    end
                end
            end
            res(j,i,k)=mol/(2*den);
        end
    end
end

for k=1:nz
    for j=1:ny
        for i=1:nx
            fprintf(fid,' %5.4f',res(j,i,k));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
        
%%%%%%%%%%%%%%%%%%%%% resolution S %%%%%%%%%%%%%%%%%%%%%%
res=ones(ny,nx,nz);

fid=fopen('res_S.dat','w');

for k=1:nz
    for i=1:nx
        for j=1:ny
            VEL(j,i,k)=vs((k-1)*ny+j,i);
            VEL_1d(j,i,k)=vs_1d((k-1)*ny+j,i);
            VEL_chk(j,i,k)=vs_checkbd((k-1)*ny+j,i);
        end
    end
end

for k=2:nz-1
    %orig=VEL(1,1,k);
    %disp(orig);
    for j=2:ny-1
        for i=2:nx-1
            mol=0.0;
            den=0.0;
           % distx=0;disty=0;
            for jj=2:ny-1
                for ii=2:nx-1
                    distx=0;disty=0;
                    if(i<ii)
                        ilow=i;iup=ii;
                    else
                        ilow=ii;iup=i;
                    end
                    if(j<jj)
                        jlow=j;jup=jj;
                    else
                        jlow=jj;jup=j;
                    end
                    for itmp=ilow:iup-1
                        distx=distx+dx(itmp);
                    end
                    for itmp=jlow:jup-1
                        disty=disty+dy(itmp);
                    end
                    dist=sqrt(distx^2+disty^2);
                    if(dist<=radius)
                        mol=mol+((VEL_chk(jj,ii,k)-VEL_1d(jj,ii,k))+(VEL(jj,ii,k)-VEL_1d(jj,ii,k)))^2;
                        den=den+(VEL_chk(jj,ii,k)-VEL_1d(jj,ii,k))^2+(VEL(jj,ii,k)-VEL_1d(jj,ii,k))^2;
                        %disp(mol);
                        %disp(den);
                    end
                end
            end
            res(j,i,k)=mol/(2*den);
        end
    end
end

for k=1:nz
    for j=1:ny
        for i=1:nx
            fprintf(fid,' %5.4f',res(j,i,k));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
            
