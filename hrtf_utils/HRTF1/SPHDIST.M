function a=sphdist(sp1,sp2)
%function a=sphdist(sp1,sp2)
%calculates the distance of two points on a sphere
%sp1=[ele1,azi1]
%sp2=[ele2,azi2]

[l c] = size(sp2);

a=zeros(1,c);

el1=sp1(1);
az1=sp1(2);
az1=az1/180*pi;
el1=el1/180*pi;

[x1,y1,z1]=sph2cart(az1,el1,1);
v1=[x1 y1 z1];

for i=1:c 
el2=sp2(1,i);
az2=sp2(2,i);

az2=az2/180*pi;
el2=el2/180*pi;

[x2,y2,z2]=sph2cart(az2,el2,1);

v2=[x2 y2 z2];

prod_v =v1*v2';

a(i)=acos(prod_v)/pi*180;
if (a(i) <1e-4) a(i)=0; end;
end;

return
end;
