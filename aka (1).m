function dydt = aka(t,n) 

N=6.022*10^23*10^-6; 
pressureCH4=0.45*101325; 
v0=1*10^-3; 
molesCH4=pressureCH4*v0/(8.314*973); 
pressureAR=0.45*101325;
molesAR=pressureAR*v0/(8.314*973);
pressureNA=0.1*101325;%pressure of sodium in pascal
molesNA=pressureNA*v0/(8.314*973);%initial moles of sodium
ni=molesCH4+molesAR+molesNA ;%total initial moles
nutotal=n(1)+n(2)+n(3)+n(4)+n(5)+n(6)+n(7)+n(8)+n(9)+n(10)+n(11)+n(12)+n(13)+n(14)+n(15)+n(16)+molesAR;

votot=nutotal*v0/ni; % total volume at any time t since pressure is constant

kf=[1.1*10^-8 2.2*10^9 1.5*10^11 N*2.2*10^-21 N*3.2*10^-22 N*10^-23 N*9.4*10^-27 N*2.3*10^-23 N*3.1*10^-21 N*6.9*10^-13 N*9.9*10^-24 9.5*10^9 8.7*10^10 6.8*10^3 8.4*10^5 7.1*10^3 5*10^9 5.8*10^8 4.4*10^-1 1.4*10^8 2.4*10^1 4.7*10^6 N*1.9*10^-10 N*2.8*10^-10];

kb=[N*4.1*10^-10 N*1.5*10^-9  N*4.4*10^-9 N*1.9*10^-10 4.1*10^9 N*1.9*10^-12 N*6.6*10^-11 3.9*10^2 4.7*10^4 N*5.7*10^-15 N*3.8*10^-19 N*3.3*10^-9 N*7.3*10^-9 N*3*10^-9  N*2.3*10^-10 N*1.4*10^-9 N*1.5*10^-9 N*3.1*10^-9 N*3.7*10^-12 N*5.1*10^-10 N*5*10^-10 N*3.9*10^-11 1.8*10^-10 2.7*10^-4];


%{
rate=d(species)/dt and denoted as ri for each species
z1=ch4 & r1
z2=ch3* & r2
z3=h* & r3
z4=Na2 & r4
z5=Na &r5
z6=Na3   & r6
z7=NaH & r7
z8=HNaCh3 & r8
z9=Na2H & r9
z10=Na2Ch3 & r10
z11=HNa2Ch3 & r11
z12=HNa3Ch3 & r12
z13=h2 & r13
z14=C2H6 & r14
z15=NaCh3 & r15
z16=Na3H & r16
%}

z1= kf(1)*(n(1)/votot) - kb(1)*(n(2)/votot)*(n(3)/votot);

z2= kf(2)*(n(4)/votot) - kb(2)*(n(5)/votot);

z3= kf(3)*(n(6)/votot)  - kb(3)*(n(4)/votot)*(n(5)/votot);

z4=  kf(4)*(n(1)/votot)*(n(5)/votot) - kb(4)*(n(7)/votot)*(n(2)/votot);

z5= kf(5)*(n(5)/votot)*(n(1)/votot) - kb(5)*(n(8)/votot);

z6= kf(6)*(n(4)/votot)*(n(1)/votot) - kb(6)*(n(9)/votot)*(n(2)/votot);

z7= kf(7)*(n(1)/votot)*(n(4)/votot) - kb(7)*(n(10)/votot)*(n(3)/votot);

z8= kf(8)*(n(1)/votot)*(n(4)/votot)  - kb(8)*(n(11)/votot);

z9= kf(9)*(n(1)/votot)*(n(6)/votot) - kb(9)*(n(12)/votot);

z10=  kf(10)*(n(1)/votot)*(n(3)/votot)  - kb(10)*(n(2)/votot)*(n(13)/votot);

z11=  kf(11)*(n(1)/votot)*(n(2)/votot) - kb(11)*(n(14)/votot)*(n(3)/votot);

z12= kf(12)*(n(8)/votot) - kb(12)*(n(15)/votot)*(n(3)/votot);

z13= kf(13)*(n(8)/votot) - kb(13)*(n(2)/votot)*(n(7)/votot);

z14=  kf(14)*(n(7)/votot) - kb(14)*(n(3)/votot)*(n(5)/votot);

z15=  kf(15)*(n(15)/votot) - kb(15)*(n(2)/votot)*(n(5)/votot);

z16=  kf(16)*(n(9)/votot) - kb(16)*(n(4)/votot)*(n(3)/votot);

z17= kf(17)*(n(9)/votot) - kb(17)*(n(5)/votot)*(n(7)/votot);

z18= kf(18)*(n(10)/votot)  - kb(18)*(n(4)/votot)*(n(2)/votot);

z19= kf(19)*(n(10)/votot)  - kb(19)*(n(5)/votot)*(n(15)/votot);

z20= kf(20)*(n(16)/votot)  - kb(20)*(n(4)/votot)*(n(7)/votot);

z21=kf(21)*(n(16)/votot)  - kb(21)*(n(6)/votot)*(n(3)/votot);

z22= kf(22)*(n(16)/votot)  - kb(22)*(n(5)/votot)*(n(9)/votot);

z23= kf(23)*(n(3)/votot)*(n(3)/votot) - kb(23)*(n(13)/votot);

z24= kf(24)*(n(2)/votot)*(n(2)/votot) - kb(24)*(n(14)/votot);

r1= - votot*(z1+z4+z5+z6+z7+z8+z9+z10+z11);

r2= votot*(z1+z4+z6+z10-z11+z13+z15+z18-2*z24);

r3= votot*(z1+z7-z10+z11+z12+z14+z16+z21-2*z23);

r4=votot*(-z2+z3-z6-z7-z8+z16+z18+z20);

r5=votot*(2*z2+z3-z4-z5+z14+z15+z17+z19+z22);

r6=votot*(-z3-z9+z21);

r7=votot*(z4+z13-z14+z17+z20);

r8=votot*(z5-z12-z13);

r9=votot*(z6-z16-z17+z22);

r10=votot*(z7-z18-z19);

r11=votot*(z8);

r12=votot*(z9);

r13=votot*(z10+z23);

r14=votot*(z11+z24);

r15=votot*(z12-z15+z19);

r16=votot*(-z20-z21-z22);

dydt=[r1; r2; r3; r4; r5; r6; r7; r8; r9; r10; r11; r12; r13; r14; r15; r16]


end
