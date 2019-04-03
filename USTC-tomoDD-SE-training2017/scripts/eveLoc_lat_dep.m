%
% plot earhtquake location in lat-dep plane.
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
for i=1:length(z)
    markersize=6.0;
    plot(y(i),z(i),'.','markersize',markersize,'color','k','linewidth',6);
    hold on
end
%axis equal;
axis ij;
axis image;
axis([Y(2),Y(ny-1),Z(2),Z(nz-4)]);
set(gca,'DataAspectRatio',[1 111 1])
title(['Relocated Earthquake Location'],'fontsize',10,'fontweight','bold');
%title(['Initial Earthquake Epicenters'],'fontsize',10,'fontweight','bold');
xlabel('Lat(degree)','fontsize',10,'fontweight','bold');
ylabel('Depth(km)','fontsize',10,'fontweight','bold');
print(h,'-depsc2','lat-dep_reloc','-r300');
