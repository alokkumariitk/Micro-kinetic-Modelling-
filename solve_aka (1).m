t=[0 7200]; %seconds

N=6.022*10^23*10^-6; 
pressureCH4=0.45*101325; 
v0=1*10^-3; 
molesCH4=pressureCH4*v0/(8.314*973); 
pressureAR=0.45*101325;
molesAR=pressureAR*v0/(8.314*973);
pressureNA=0.1*101325;%pressure of sodium in pascal
molesNA=pressureNA*v0/(8.314*973);%initial moles of sodium
ni=molesCH4+molesAR+molesNA %total initial moles
nutotal=n(1)+n(2)+n(3)+n(4)+n(5)+n(6)+n(7)+n(8)+n(9)+n(10)+n(11)+n(12)+n(13)+n(14)+n(15)+n(16)+molesAR;

votot=nutotal*v0/ni % total volume at any time t since pressure is constant

j=[nch4 0 0 0 nna 0 0 0 0 0 0 0 0 0 0 0]
[t,n]=ode15s(@aka,t,j)
nutotal=n(:,1)+n(:,2)+n(:,3)+n(:,4)+n(:,5)+n(:,6)+n(:,7)+n(:,8)+n(:,9)+n(:,10)+n(:,11)+n(:,12)+n(:,13)+n(:,14)+n(:,15)+n(:,16)+nar;
vt=nutotal*v0/ni;
kf=[1.1*10^-8 2.2*10^9 1.5*10^11 N*2.2*10^-21 N*3.2*10^-22 N*10^-23 N*9.4*10^-27 N*2.3*10^-23 N*3.1*10^-21 N*6.9*10^-13 N*9.9*10^-24 9.5*10^9 8.7*10^10 6.8*10^3 8.4*10^5 7.1*10^3 5*10^9 5.8*10^8 4.4*10^-1 1.4*10^8 2.4*10^1 4.7*10^6 N*1.9*10^-10 N*2.8*10^-10];

kb=[N*4.1*10^-10 N*1.5*10^-9  N*4.4*10^-9 N*1.9*10^-10 4.1*10^9 N*1.9*10^-12 N*6.6*10^-11 3.9*10^2 4.7*10^4 N*5.7*10^-15 N*3.8*10^-19 N*3.3*10^-9 N*7.3*10^-9 N*3*10^-9  N*2.3*10^-10 N*1.4*10^-9 N*1.5*10^-9 N*3.1*10^-9 N*3.7*10^-12 N*5.1*10^-10 N*5*10^-10 N*3.9*10^-11 1.8*10^-10 2.7*10^-4];
z1= kf(1)*(n(:,1)/vt) - kb(1)*(n(:,2)/vt).*(n(:,3)/vt);

z2= kf(2)*(n(:,4)/vt) - kb(2)*(n(:,5)/vt);

z3= kf(3)*(n(:,6)/vt)  - kb(3)*(n(:,4)/vt).*(n(:,5)/vt);

z4=  kf(4)*(n(:,1)/vt).*(n(:,5)/vt) - kb(4)*(n(:,7)/vt).*(n(:,2)/vt);

z5= kf(5)*(n(:,5)/vt).*(n(:,1)/vt) - kb(5)*(n(:,8)/vt);

z6= kf(6)*(n(:,4)/vt).*(n(:,1)/vt) - kb(6)*(n(:,9)/vt).*(n(:,2)/vt);

z7= kf(7)*(n(:,1)/vt).*(n(:,4)/vt) - kb(7)*(n(:,10)/vt).*(n(:,3)/vt);

z8= kf(8)*(n(:,1)/vt).*(n(:,4)/vt)  - kb(8)*(n(:,11)/vt);

z9= kf(9)*(n(:,1)/vt).*(n(:,6)/vt) - kb(9)*(n(:,12)/vt);

z10=  kf(10)*(n(:,1)/vt).*(n(:,3)/vt)  - kb(10)*(n(:,2)/vt).*(n(:,13)/vt);

z11=  kf(11)*(n(:,1)/vt).*(n(:,2)/vt) - kb(11)*(n(:,14)/vt).*(n(:,3)/vt);

z12= kf(12)*(n(:,8)/vt) - kb(12)*(n(:,15)/vt).*(n(:,3)/vt);

z13= kf(13)*(n(:,8)/vt) - kb(13)*(n(:,2)/vt).*(n(:,7)/vt);

z14=  kf(14)*(n(:,7)/vt) - kb(14)*(n(:,3)/vt).*(n(:,5)/vt);

z15=  kf(15)*(n(:,15)/vt) - kb(15)*(n(:,2)/vt).*(n(:,5)/vt);

z16=  kf(16)*(n(:,9)/vt) - kb(16)*(n(:,4)/vt).*(n(:,3)/vt);

z17= kf(17)*(n(:,9)/vt) - kb(17)*(n(:,5)/vt).*(n(:,7)/vt);

z18= kf(18)*(n(:,10)/vt)  - kb(18)*(n(:,4)/vt).*(n(:,2)/vt);

z19= kf(19)*(n(:,10)/vt)  - kb(19)*(n(:,5)/vt).*(n(:,15)/vt);

z20= kf(20)*(n(:,16)/vt)  - kb(20)*(n(:,4)/vt).*(n(:,7)/vt);

z21=kf(21)*(n(:,16)/vt)  - kb(21)*(n(:,6)/vt).*(n(:,3)/vt);

z22= kf(22)*(n(:,16)/vt)  - kb(22)*(n(:,5)/vt).*(n(:,9)/vt);

z23= kf(23)*(n(:,3)/vt).*(n(:,3)/vt) - kb(23)*(n(:,13)/vt);

z24= kf(24)*(n(:,2)/vt).*(n(:,2)/vt) - kb(24)*(n(:,14)/vt);

r1= - votot.*(z1(:,:)+z4(:,:)+z5(:,:)+z6(:,:)+z7(:,:)+z8(:,:)+z9(:,:)+z10(:,:)+z11(:,:))
plot(t,n(:,1))
figure

x1= (kf(1).*n(:,1)-kb(1).*(n(:,2)./(vt)).*(n(:,3)./(vt))).*(vt./r1);
x4=(kf(4).*(n(:,1)).*(n(:,5)./vt)-kb(4).*(n(:,7)./(vt)).*(n(:,2)./(vt))).*(vt./r1);
x5=(kf(5).*(n(:,1)./(vt)).*(n(:,5)./vt)-kb(5).*(n(:,8)./(vt))).*(vt./r1);
x6=(kf(6).*(n(:,1)./(vt)).*(n(:,4)./vt)-kb(6).*(n(:,7)./(vt)).*(n(:,2)./(vt))).*(vt./r1);
x7=(kf(7).*(n(:,1)./(vt)).*(n(:,4)./vt)-kb(7).*(n(:,10)./(vt)).*(n(:,3)./(vt))).*(vt./r1);
x8=(kf(8).*(n(:,1)./(vt)).*(n(:,4)./vt)-kb(8).*(n(:,11)./(vt))).*(vt./r1);
x9=(kf(9).*(n(:,1)./(vt)).*(n(:,6)./vt)-kb(9).*(n(:,12)./(vt))).*(vt./r1);
x10=(kf(10).*(n(:,1)./(vt)).*(n(:,3)./vt)-kb(10).*(n(:,2)./(vt)).*(n(:,13)./(vt))).*(vt./r1);
x11=(kf(11).*(n(:,1)./(vt)).*(n(:,2)./vt)-kb(11).*(n(:,14)./(vt)).*(n(:,3)./(vt))).*(vt./r1);

x2=0;
x3=0;
x12=0;
x13=0;
x14=0;
x15=0;
x16=0;
x17=0;
x18=0;
x19=0;
x20=0;
x21=0;
x22=0;
x23=0;
x24=0;

t1=0:1:14;
plot(t1,x4(1:15),t1,x5(1:15),t1,x6(1:15),t1,x8(1:15),t1,x9(1:15),t1,x7(1:15),t1,x1(1:15),t1,x10(1:15),t1,x11(1:15))