%
% plot earhtquake location in lat-lon plane.
%
clear;

%%%%%%%%%%%%%%%%%% read grid %%%%%%%%%%%%%%%%%
fid=fopen('../MOD','r');
gridnumber=fscanf(fid,'%f',[1,4]);
nx=gridnumber(2);ny=gridnumber(3);nz=gridnumber(4);
X=fscanf(fid,'%f',[1,nx]);
Y=fscanf(fid,'%f',[1,ny]);
Z=fscanf(fid,'%f',[1,nz]);

%%%%%%%%%%%%%%%%%% read earthquake location %%%%%%%%%%%%%
file_loc='tomoDD.reloc';
mdat=load(file_loc); 
cusp = mdat(:,1); mag= mdat(:,17); lon=mdat(:,3); lat=mdat(:,2);
%x = mdat(:,5)/1000; y = mdat(:,6)/1000; z = mdat(:,4);
x = mdat(:,3); y = mdat(:,2); z = mdat(:,4);
ex = mdat(:,8)/1000; ey = mdat(:,9)/1000; ez = mdat(:,10)/1000;
yr = mdat(:,11);mo = mdat(:,12); dy = mdat(:,13);

%%%%%%%%%%%%%%%%%% plot earthquake location %%%%%%%%%%%%%%%%
h=figure;
%subplot(2,2,1)
for i=1:length(z)
    markersize = 6;
    plot(x(i),y(i),'.','markersize',markersize,'color','k','linewidth',6);
    hold on
end
% plot station
fid=fopen('../Input_Files/station.dat','r');
while ~feof(fid)
    [sta]=fscanf(fid,'%s',1);[sta_loc]=fscanf(fid,'%f',[1,3]);
    plot(sta_loc(2),sta_loc(1),'r^');
    if sta_loc(1)>Y(2) && sta_loc(1)<Y(ny-1) && sta_loc(2)>X(2) && sta_loc(2)<X(nx-1)
        %text(sta_loc(2),sta_loc(1)+0.1,sta,'fontsize',8,'color','r');
    end
    fscanf(fid,'\n');
end
axis([X(2),X(nx-1),Y(2),Y(ny-1)]);
%set(gca,'DataAspectRatio',[5 1 1])
title(['Relocated Earthquake Epicenters'],'fontsize',10,'fontweight','bold');
%title(['Initial Earthquake Epicenters'],'fontsize',10,'fontweight','bold');
xlabel('Lon(degree)','fontsize',10,'fontweight','bold');
ylabel('Lat(degree)','fontsize',10,'fontweight','bold');
print(h,'-depsc2','lonlat_reloc','-r300');
