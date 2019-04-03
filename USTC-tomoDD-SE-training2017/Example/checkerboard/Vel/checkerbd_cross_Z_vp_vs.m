%
% plot Z (depth) cross sections of vp and vs.
% Hao Guo.
%
clc
clear
%close all

%%%%%%%%%%%%%%%%%%%% Parameter %%%%%%%%%%%%%%%%%%%%%%%%
pinter = 0; % 1:interpolate the velocity; 0: dont interpolate.
if pinter==1
    inter1 = 0.05; % interpolation in Y direction
    inter2 = 0.05;    % interpolation in Z direction
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
fid=fopen('MOD','r');
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

[x1,y1]=meshgrid(X(2:nx-1),Y(2:ny-1));

%%%%%%%%%%%%%  Draw the cross-sections at Y=Y(j) km %%%%%%%%%%%%%%%%%%%
min_X = X(2); max_X = X(nx-1);  % min and max limit of X (lon)
min_Y = Y(2); max_Y = Y(ny-1);      % min and max limit of Z

for k=2:nz-1
    
    %%%%%%%%%%%%%%%%% subplot Vp %%%%%%%%%%%%%%%%%
    h=figure;
    subplot(2,1,1);
    for j=2:ny-1
        for i=2:nx-1
            crossh(j-1,i-1)=VEL_P(i,j,k);
            dwsh(j-1,i-1)=DWS_P(i,j,k);
            resh(j-1,i-1)=RES_P(i,j,k);
        end
    end
    
    caxis([-ano_Vp,ano_Vp]);
    caxis manual;
    colormap(ColorNJet)
    hold on
    
    if pinter==1
        [xi,yi] = meshgrid(min_X:inter1:max_X, min_Y:inter2:max_Y);
        cross_inter = interp2(x1,y1,crossh,xi,yi,'linear');
        
        if pres==1
            dws_inter_P = interp2(x1,y1,dwsh,xi,yi,'linear');
            res_inter_P = interp2(x1,y1,resh,xi,yi,'linear');
            
            s=size(cross_inter);
            for j=1:s(1)
                for i=1:s(2)
                    if dws_inter_P(j,i)<(mean(mean(dws_inter_P)))*dwsthres || res_inter_P(j,i)<resthres
                        cross_inter(j,i)=nan;
                    end
                end
            end
        end
        
        pcolor(xi,yi,cross_inter);
        shading flat;
    else
        pcolor(x1,y1,crossh);
    end
    
    for i=1:length(z)
        if abs(y(i))>90
            y(i) = y(i)-y(i)/abs(y(i))*180;
        end
        if abs(x(i))>180
            x(i) = x(i)-x(i)/abs(x(i))*360;
        end
        if k>=2 && k<=nz-1 && -z(i)>(Z(k)+Z(k-1))/2 && -z(i)<=(Z(k)+Z(k+1))/2
            plot(x(i),y(i),'.','markersize',makersize,'color','k');
        end
    end
    
    title(['Vp Z=',num2str(Z(k)),' km'],'fontsize',12,'fontweight','bold');%%%%%%%%%%%%%%%%%%%%%%%
    xlabel('Lon(degree)','fontsize',12,'fontweight','bold');
    ylabel('lat(degree)','fontsize',12,'fontweight','bold');
    axis image;
    axis ij;
    axis([min_X,max_X,min_Y,max_Y]);
    my_handle=colorbar;
    %set(get(my_handle,'Title'),'string','km/s');
    
    
    %%%%%%%%%%%%%%%%% subplot Vs %%%%%%%%%%%%%%%%%
    subplot(2,1,2);
    
    for j=2:ny-1
        for i=2:nx-1
            crossh(j-1,i-1)=VEL_S(i,j,k);
            dwsh(j-1,i-1)=DWS_S(i,j,k);
            resh(j-1,i-1)=RES_S(i,j,k);
        end
    end
    
    caxis([-ano_Vp,ano_Vp]);
    caxis manual;
    colormap(ColorNJet)
    hold on
    if pinter==1
        [xi,yi] = meshgrid(min_X:inter1:max_X, min_Y:inter2:max_Y);
        cross_inter = interp2(x1,y1,crossh,xi,yi,'linear');
        
        if pres==1
            dws_inter_S = interp2(x1,y1,dwsh,xi,yi,'linear');
            res_inter_S = interp2(x1,y1,resh,xi,yi,'linear');
            
            s=size(cross_inter);
            for j=1:s(1)
                for i=1:s(2)
                    if dws_inter_S(j,i)<(mean(mean(dws_inter_S)))*dwsthres || res_inter_S(j,i)<resthres
                        cross_inter(j,i)=nan;
                    end
                end
            end
        end
        
        pcolor(xi,yi,cross_inter);
        shading flat;
    else
        pcolor(x1,y1,crossh);
    end
    for i=1:length(z)
        if abs(y(i))>90
            y(i) = y(i)-y(i)/abs(y(i))*180;
        end
        if abs(x(i))>180
            x(i) = x(i)-x(i)/abs(x(i))*360;
        end
        if k>=2 && k<=nz-1 && -z(i)>(Z(k)+Z(k-1))/2 && -z(i)<=(Z(k)+Z(k+1))/2
            plot(x(i),y(i),'.','markersize',makersize,'color','k');
        end
    end
    
    title(['Vs Z=',num2str(Z(k)),' km'],'fontsize',12,'fontweight','bold');%%%%%%%%%%%%%%%%%%%%%%%
    xlabel('Lon(degree)','fontsize',12,'fontweight','bold');
    ylabel('lat(degree)','fontsize',12,'fontweight','bold');
    axis image;
    axis ij;
    axis([min_X,max_X,min_Y,max_Y]);
    my_handle=colorbar;
    %set(get(my_handle,'Title'),'string','km/s');
    
    %print(h,'-dpng',strcat('inter',num2str(pinter),'_vp&vs_Z',num2str(Z(k)),'_lon',num2str(min_X),'-',num2str(max_X),'_lat',num2str(min_Y),'-',num2str(max_Y)),'-r300');
end
