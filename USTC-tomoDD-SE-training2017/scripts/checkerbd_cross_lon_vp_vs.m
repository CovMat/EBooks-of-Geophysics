%
% plot X (longitude) cross sections of vp and vs.
% Hao Guo.
%
clc
clear
%close all

%%%%%%%%%%%%%%%%%%%% Parameter %%%%%%%%%%%%%%%%%%%%%%%%
pinter = 1; % 1:interpolate the velocity; 0: dont interpolate.
if pinter==1
    inter1 = 0.01; % interpolation in Y direction
    inter2 = 1;    % interpolation in Z direction
    pres = 0;      % 1:use prepared resolution files to cut off regions with bad resolution; 0: dont use.
    if pres==1
        dwsthres = 0.02; % dws threshold value
        resthres = 0.8;  % resolution threshold value
    end
end

makersize = 5; % earthquake markersize

ano_Vp = 0.05;
ano_Vs = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% read relocation file %%%%%%%%%%%%%%%%%%%
file_loc='tomoDD.reloc';
mdat=load(file_loc);
cusp = mdat(:,1); mag= mdat(:,17); lon=mdat(:,3); lat=mdat(:,2);
%x = mdat(:,5)/1000; y = mdat(:,6)/1000; z = -mdat(:,4);
x = mdat(:,3); y = mdat(:,2); z = -mdat(:,4);
ex = mdat(:,8)/1000; ey = mdat(:,9)/1000; ez = mdat(:,10)/1000;
yr = mdat(:,11);mo = mdat(:,12); dy = mdat(:,13);

%%%%%%%%%%%%%%%%%%%% real initial velocity model %%%%%%%%%%%%%%%%%%
fid=fopen('../MOD','r');
gridnumber=fscanf(fid,'%f',[1,4]);
nx=gridnumber(2);ny=gridnumber(3);nz=gridnumber(4);
X=fscanf(fid,'%f',[1,nx]);
Y=fscanf(fid,'%f',[1,ny]);
Z=fscanf(fid,'%f',[1,nz]);
for k=1:nz
    for j=1:ny
        for i=1:nx
            vp_ini((k-1)*ny+j,i)=fscanf(fid,'%f',1);
        end
    end
end
for k=1:nz
    for j=1:ny
        for i=1:nx
            vpvs_ini((k-1)*ny+j,i)=fscanf(fid,'%f',1);
        end
    end
end
vs_ini = vp_ini./vpvs_ini;

%%%%%%%%%%%%%%%%%%%%%%%%%% read inverted velocity data %%%%%%%%%%%%%%%%%
vp=load('Vp_model.dat');
vs=load('Vs_model.dat');

%%%%%%%%%%%%%%%%%%%%%%%%%% read resolution files %%%%%%%%%%%%%%%%%%%%%%%%%%
if pinter==1 && pres==1
    dws_P=load('DWS_P');
    res_P=load('res_P.dat');
    dws_S=load('DWS_S');
    res_S=load('res_S.dat');
else
    dws_P=ones(ny*nz,nx);
    res_P=ones(ny*nz,nx);
    dws_S=ones(ny*nz,nx);
    res_S=ones(ny*nz,nx);
end
%%%%%%%%%%%%%%%%%%%%%%%% 2D model to 3D model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:nz
    for j=1:ny
        for i=1:nx
            % Vp
            VEL_P(i,j,k)=(vp((k-1)*ny+j,i)-vp_ini((k-1)*ny+j,i))/vp_ini((k-1)*ny+j,i);
            DWS_P(i,j,k)=dws_P((k-1)*ny+j,i);
            RES_P(i,j,k)=res_P((k-1)*ny+j,i);
            % Vs
            VEL_S(i,j,k)=(vs((k-1)*ny+j,i)-vs_ini((k-1)*ny+j,i))/vs_ini((k-1)*ny+j,i);
            DWS_S(i,j,k)=dws_S((k-1)*ny+j,i);
            RES_S(i,j,k)=res_S((k-1)*ny+j,i);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% set colormap %%%%%%%%%%%%%%%%%%%%%%%
if pinter==1
    ColorJet=colormap('Jet');
elseif pinter==0
    ColorJet=colormap('gray');
end
ColorNJet=flipud(ColorJet);

[y1,z1]=meshgrid(Y(2:ny-1),Z(2:nz-1));

%%%%%%%%%%%%%  Draw the cross-sections at X=X(i) km %%%%%%%%%%%%%%%%%%%
min_Y = Y(2); max_Y = Y(ny-1);  % min and max limit of Y (lat)
min_Z = Z(2); max_Z =Z(nz-1);      % min and max limit of Z

for i=2:nx-1
    
    %%%%%%%%%%%%%%%%% subplot Vp %%%%%%%%%%%%%%%%%
    h=figure;
    subplot(2,1,1);
    for k=2:nz-1
        for j=2:ny-1
            crossh(k-1,j-1)=VEL_P(i,j,k);
            dwsh(k-1,j-1)=DWS_P(i,j,k);
            resh(k-1,j-1)=RES_P(i,j,k);
        end
    end
    
    caxis([-ano_Vp,ano_Vp]);
    caxis manual;
    colormap(ColorNJet)
    hold on
    
    if pinter==1
        
        [yi,zi] = meshgrid(min_Y:inter1:max_Y, min_Z:inter2:max_Z);
        cross_inter = interp2(y1,z1,crossh,yi,zi,'linear');
        
        if pres==1
            dws_inter_P = interp2(y1,z1,dwsh,yi,zi,'linear');
            res_inter_P = interp2(y1,z1,resh,yi,zi,'linear');
            
            s=size(cross_inter);
            for k=1:s(1)
                for j=1:s(2)
                    if dws_inter_P(k,j)<(mean(mean(dws_inter_P)))*dwsthres || res_inter_P(k,j)<resthres
                        cross_inter(k,j)=nan;
                    end
                end
            end
        end
        
        pcolor(yi,zi,cross_inter);
        shading flat;
    else
        pcolor(y1,z1,crossh);
    end
    
    for j=1:length(z)
        if abs(y(j))>90
            y(j) = y(j)-y(j)/abs(y(j))*180;
        end
        if abs(x(j))>180
            x(j) = x(j)-x(j)/abs(x(j))*360;
        end
        if i>=2 && x(j)>(X(i)+X(i-1))/2 && x(j)<=(X(i)+X(i+1))/2
            plot(y(j),-z(j),'.','markersize',makersize,'color','k');
        end
    end
    
    title(['Vp Lon=',num2str(X(i)),' degree'],'fontsize',12,'fontweight','bold');%%%%%%%%%%%%%%%%%%%%%%%
    xlabel('Lat(degree)','fontsize',12,'fontweight','bold');
    ylabel('Z(km)','fontsize',12,'fontweight','bold');
    axis image;
    axis ij;
    axis([min_Y,max_Y,min_Z,max_Z]);
    my_handle=colorbar;
    %set(get(my_handle,'Title'),'string','km/s');
    set(gca,'DataAspectRatio',[1 111 1]);
    
    
    %%%%%%%%%%%%%%%%% subplot Vs %%%%%%%%%%%%%%%%%
    
    subplot(2,1,2);
    for k=2:nz-1
        for j=2:ny-1
            crossh(k-1,j-1)=VEL_S(i,j,k);
            dwsh(k-1,j-1)=DWS_S(i,j,k);
            resh(k-1,j-1)=RES_S(i,j,k);
        end
    end
    
    caxis([-ano_Vp,ano_Vp]);
    caxis manual;
    colormap(ColorNJet)
    hold on
    
    if pinter==1
        
        [yi,zi] = meshgrid(min_Y:inter1:max_Y, min_Z:inter2:max_Z);
        cross_inter = interp2(y1,z1,crossh,yi,zi,'linear');
        
        if pres==1
            dws_inter_S = interp2(y1,z1,dwsh,yi,zi,'linear');
            res_inter_S = interp2(y1,z1,resh,yi,zi,'linear');
            
            s=size(cross_inter);
            for k=1:s(1)
                for j=1:s(2)
                    if dws_inter_S(k,j)<(mean(mean(dws_inter_S)))*dwsthres || res_inter_S(k,j)<resthres
                        cross_inter(k,j)=nan;
                    end
                end
            end
        end
        
        pcolor(yi,zi,cross_inter);
        shading flat;
    else
        pcolor(y1,z1,crossh);
    end
    for j=1:length(z)
        if abs(y(j))>90
            y(j) = y(j)-y(j)/abs(y(j))*180;
        end
        if abs(x(j))>180
            x(j) = x(j)-x(j)/abs(x(j))*360;
        end
        if i>=2 && x(j)>(X(i)+X(i-1))/2 && x(j)<=(X(i)+X(i+1))/2
            plot(y(j),-z(j),'.','markersize',makersize,'color','k');
        end
    end
    
    title(['Vs Lon=',num2str(X(i)),' degree'],'fontsize',12,'fontweight','bold');%%%%%%%%%%%%%%%%%%%%%%%
    xlabel('Lat(degree)','fontsize',12,'fontweight','bold');
    ylabel('Z(km)','fontsize',12,'fontweight','bold');
    axis image;
    axis ij;
    axis([min_Y,max_Y,min_Z,max_Z]);
    my_handle=colorbar;
    %set(get(my_handle,'Title'),'string','km/s');
    set(gca,'DataAspectRatio',[1 111 1]);
    
    print(h,'-dpng',strcat('inter',num2str(pinter),'_vp&vs_lon',num2str(X(i)),'_lat',num2str(min_Y),'-',num2str(max_Y),'_Z',num2str(min_Z),'-',num2str(max_Z)),'-r300');
end
