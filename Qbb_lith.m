function Qbb_lith = Qbb_lith(r0,nspoke,lambda,n)
%input parameters are: 
%r0 as in Fig. 2 [microns], this is the radius of the coccosphere which relates to the size of the coccoliths as given below 
%nspoke=number of spokes in the coccolith [dim], 
%lambda=light wavelength in water [microns], 
%n=real part of the refractive index of the coccolith material

r1=0.18*r0;%radius of central cylinder joining the proximal and distal sheets
r2=0.46*r0;%radius of the proximal sheet 
r3=0.54*r0;%radius of the distal sheet 
rm=0.50*r0;%mean of the distal and proximal sheet lateral dimensions which we will use for the equivalent disk, =(r3+r2)/2
dc=0.07*r0;%thickness of proximal and distal sheets
dh=0.18*r0;%center to center spacing of the proximal and distal sheets
rtot=1.25*r0;%maximum radius of the spherical approximation to a coccolithophore
tt=0.164*r0;%equivalent disk thickness; a disk of area Pi rm^2 with a thickness tt has the same volume as a Zhai coccolith

x0=2*pi*r0/lambda;
rmaxspoke=(0.54-0.07/2)*r0;
xmaxspoke=(0.54-0.07/2)*x0;
rminspoke=(0.18+0.07/2)*r0;
xminspoke=(0.18+0.07/2)*x0;
x0min=(nspoke/2)/(0.54-0.07/2);
x0max=(nspoke/2)/(0.18+0.07/2);
rlith=(0.54+0.07/2)*r0;
xlith=(0.54+0.07/2)*x0;
%F_rg is the maximum fraction of the projected area of the particle that can become rough
F_rg=(xmaxspoke^2-xminspoke^2)/xlith^2; % Eq. (9)
%R is a mixing coefficient
R=(1-((nspoke/2)/(xmaxspoke))^2)/(1-(xminspoke/xmaxspoke)^2);% Eq. (8) 
if isreal(n)
    %The factor 4 below corresponds to the number of scattering surfaces in a coccolith (see Fig. A1)
        Qbb_lith_trans=4*omega_sm(n)*(1-F_rg*R)+4*omega_rg(n)*F_rg*R; % Eq. (3)
        Qbb_lith_R0=4*omega_sm(n);%R=0;
        Qbb_lith_R1=4*omega_sm(n)*(1-F_rg)+4*omega_rg(n)*F_rg;%R=1;
        if x0<x0min
            Qbb_lith=Qbb_lith_R0;
        elseif x0>x0max
            Qbb_lith=Qbb_lith_R1;
        else
            Qbb_lith=Qbb_lith_trans;
        end
else
    k=imag(n);
    m=real(n);
    neff=1+sqrt((m-1)^2+k^2);
    tt=0.164*r0;%equivalent disk thickness
    %pure water in the gap between the distal and proximal sheets
    %Pope and Fry pure water absorption values
    wl_water=[380	382.5	385	387.5	390	392.5	395	397.5	400	402.5	405	407.5	410	412.5	415	417.5	420	422.5	425	427.5	430	432.5	435	437.5	440	442.5	445	447.5	450	452.5	455	457.5	460	462.5	465	467.5	470	472.5	475	477.5	480	482.5	485	487.5	490	492.5	495	497.5	500	502.5	505	507.5	510	512.5	515	517.5	520	522.5	525	527.5	530	532.5	535	537.5	540	542.5	545	547.5	550	552.5	555	557.5	560	562.5	565	567.5	570	572.5	575	577.5	580	582.5	585	587.5	590	592.5	595	597.5	600	602.5	605	607.5	610	612.5	615	617.5	620	622.5	625	627.5	630	632.5	635	637.5	640	642.5	645	647.5	650	652.5	655	657.5	660	662.5	665	667.5	670	672.5	675	677.5	680	682.5	685	687.5	690	692.5	695	697.5	700	702.5	705	707.5	710	712.5	715	717.5	720	722.5	725	727.5	730  732.5  735.0  737.5  740.0  742.5  745.0  747.5  750.0  752.5  755.0  757.5  760.0  762.5  765.0  767.5  770.0  772.5  775.0  777.5 780.0  782.5  785.0  787.5  790.0  792.5  795.0  797.5  800.0];
    a_water=[0.01137 0.010044 0.00941 0.00917 0.00851 0.00829 0.00813 0.00775 0.00663 0.00579 0.0053 0.00503 0.00473 0.00452 0.00444 0.00442 0.00454 0.00474 0.00478 0.00482 0.00495 0.00504 0.0053	0.0058 0.00635	0.00696	0.00751	0.0083 0.00922	0.00969 0.00962 0.00957 0.00979 0.01005 0.01011 0.0102 0.0106 0.0109 0.0114 0.0121 0.0127 0.0131 0.0136 0.0144 0.015 0.0162 0.0173 0.0191 0.0204 0.0228 0.0256 0.028 0.0325 0.0372 0.0396 0.0399 0.0409 0.0416 0.0417 0.0428 0.0434 0.0447 0.0452 0.0466 0.0474 0.0489 0.0511 0.0537 0.0565 0.0593 0.0596 0.0606 0.0619 0.064 0.0642 0.0672 0.0695 0.0733 0.0772 0.0836 0.0896 0.0989 0.11 0.122 0.1351 0.1516 0.1672 0.1925 0.2224 0.247 0.2577 0.2629 0.2644 0.2665 0.2678 0.2707 0.2755 0.281 0.2834 0.2904 0.2916 0.2995 0.3012 0.3077 0.3108 0.322 0.325 0.335 0.34 0.358 0.371 0.393 0.41 0.424 0.429 0.436 0.439 0.448 0.448 0.461 0.465 0.478 0.486 0.502 0.516 0.538 0.559 0.592 0.624	0.663 0.704 0.756 0.827 0.914 1.007 1.119 1.231 1.356 1.489 1.678 1.7845 1.9333 2.0822 2.2311 2.3800 2.4025 2.4250 2.4475 2.4700 2.4900 2.5100 2.5300 2.5500 2.5400 2.5300 2.5200 2.5100 2.4725 2.4350 2.3975 2.3600 2.3100 2.2600 2.2100 2.1600 2.1375 2.1150 2.0925 2.0700];
    k_water=a_water.*wl_water*10^-9/(4*pi);%wl_water is in nanometers
    kg=interp1(wl_water,k_water,lambda*1000);%lambda is in microns
    if lambda<0.38 
        kg=10^-9;
    end
    omega_sm_a=omega_sm(neff)+omega_sm(neff)*(1-omega_sm(neff))*(1-Qabs(k,tt/2,lambda))^2+...
        omega_sm(neff)*(1-omega_sm(neff))^2*(1-Qabs(k,tt/2,lambda))^2*(1-Qabs(kg,dh,lambda))^2+...
        omega_sm(neff)*(1-omega_sm(neff))^3*(1-Qabs(k,tt/2,lambda))^4*(1-Qabs(kg,dh,lambda))^2;
    omega_rg_a=omega_rg(neff)+omega_rg(neff)*(1-omega_rg(neff))*(1-Qabs(k,tt/2,lambda))^2+...
        omega_rg(neff)*(1-omega_rg(neff))^2*(1-Qabs(k,tt/2,lambda))^2*(1-Qabs(kg,dh,lambda))^2+...
        omega_rg(neff)*(1-omega_rg(neff))^3*(1-Qabs(k,tt/2,lambda))^4*(1-Qabs(kg,dh,lambda))^2;
    Qbb_lith_trans=omega_sm_a*(1-F_rg*R)+omega_rg_a*F_rg*R; % Eq. (3)
    Qbb_lith_R0=omega_sm_a;%R=0;
    Qbb_lith_R1=omega_sm_a*(1-F_rg)+omega_rg_a*F_rg;%R=1;
    if x0<x0min
        Qbb_lith=Qbb_lith_R0;
    elseif x0>x0max
        Qbb_lith=Qbb_lith_R1;
    else
        Qbb_lith=Qbb_lith_trans;
    end
end
end
%omega_sm is the specular reflection component of the total reflection coefficient coming from the smooth part of the particle surface 
function omega_sm = omega_sm(n)
    omega_sm=(3*n^4-16*n^3+12*n^2-1+2*(2*n^2-1)^(3/2))/(6*(n^2-1)^2)*(1+3-log(16)+(37/40)*((n-1)/(n+1)))/2;
end 
%omega_rg is the diffuse reflection component originating from the rough part of the particle surface
function omega_rg = omega_rg(n)
    omega_par=1/((n^2+1)^3*(n^2-1)^2)*((n^4-1)*(n^6-4*n^5-7*n^4+4*n^3-n^2-1)+2*n^2*((n^2-1)^4*log((n-1)/(n+1))+8*n^2*(n^4+1)*log(n)));
    omega_perp=(3*n+1)*(n-1)/(3*(n+1)^2);
    omega_rg=5/6*(omega_par+omega_perp)/2;
end
function Qabs = Qabs(k,t,lambda)
    % this is an approximation to expression A4 in the paper that works if you don't have access to the exponential integral function of order 3, expint(3,x)
    Qabs=1-exp(-3*k*2*pi*t/lambda);
end