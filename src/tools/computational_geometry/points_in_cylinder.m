function I=points_in_cylinder(vc1,vc2,d,vp)
%move first coordinate of cylinder to origin
vc2=vc2-vc1;
vp=vp-vc1;
vc1=vc1-vc1;
%align the cylinder to x-axis or y-axis
x=[1 ; 0 ; 0];
if norm(cross(vc2,x))>1e-6
R=Rotate(vc2,x);
vc2=R*vc2;
vp=R*vp;
I=vp(1,:)>=vc1(1,:) & vp(1,:)<=vc2(1,:) & (vp(2,:).^2+vp(3,:).^2)<=d^2;
else
y=[0 ; 1 ; 0];
R=Rotate(vc2,y);
vc2=R*vc2;
vp=R*vp;
I=vp(2,:)>=vc1(2,:) & vp(2,:)<=vc2(2,:) & (vp(1,:).^2+vp(3,:).^2)<=d^2;
end
end

function R=Rotate(a,b)
a=a/norm(a);
b=b/norm(b);
v=cross(a,b);
s=norm(v);
c=dot(a,b);
vx=[0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R=eye(3)+vx+vx*vx*(1-c)/s^2;
end