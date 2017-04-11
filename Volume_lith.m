function vcoccolithfilled = Volume_lith(r0) 
%Zhai model dimensional relationships
dc=0.07*r0;
dh=0.18*r0;
r2=0.46*r0;
r1=0.18*r0;
r3=0.54*r0;
rm=0.50*r0;

%spherical cap volume formula, sphere radius=r, cap lateral dimension=y, cap thickness=tcap
vcap = @(r,y,tcap) 2*pi*r^2*(1-sqrt(1-(y^2/r^2)))*tcap;
%tube volume formula, outer radius=router,inner radius=rinner,height=height
vtube = @(router,rinner,height) pi*(router^2-rinner^2)*height;
%volume of filled cap ring,rsouter= outer radius of cap ring ares, rsinner= inner radius of cap ring area
vring = @(r,rsouter,rsinner,tcap) vcap(r,rsouter,tcap)-vcap(r,rsinner,tcap);
%volume of slit opening area, rsouter= outer radius of slitted ares, rsinner= inner radius of slitted area *)
vslits = @(r,rsouter,rsinner,tcap) 1/2*(vring(r,rsouter,rsinner,tcap)); 

%volume of proximal sheet of coccolith = volume of outer ring+volume of slit area+volume of inner ring
vproximal=vring(r0,r2+dc/2,r2-dc/2,dc)+vslits(r0,r2-dc/2,r1+dc/2,dc)+vring(r0,r1+dc/2,r1-dc/2,dc);
vproximalfilled=vring(r0,r2+dc/2,r2-dc/2,dc)+vring(r0,r2-dc/2,r1+dc/2,dc)+vring(r0,r1+dc/2,r1-dc/2,dc);
%volume of distal sheet of coccolith = volume of outer ring+volume od slit area+volume of inner ring
rd=r0+dh;
vdistal=vring(rd,r3+dc/2,r3-dc/2,dc)+vslits(rd,r3-dc/2,r1+dc/2,dc)+vring(rd,r1+dc/2,r1-dc/2,dc);
vdistalfilled=vring(rd,r3+dc/2,r3-dc/2,dc)+vring(rd,r3-dc/2,r1+dc/2,dc)+vring(rd,r1+dc/2,r1-dc/2,dc);
%volume of joining tube=volume of central tube-volume of inner ring of distal and proximal sheets
vcolumn=vtube(r1+dc/2,r1-dc/2,dh-dc);
%total coccolith model volume
vcoccolith=vproximal+vdistal+vcolumn;
%total coccolith model volume with slits filled
vcoccolithfilled=vproximalfilled+vdistalfilled+vcolumn;
end
