
% mP=3.242e-4;
mP1=3.141e-4;mP2=1.57e-4;

% nQ=8.83e-4;
nQ1=8.33e-4;nQ2=5.0e-4;

wc=500;
Lf=1.35e-3; Cf=50e-6; Rf=0.1; rf=0.1;
Kpv=0.123;
Kiv=1.29;
Kpi=2.025;
Kii=150;
F=0.75;
Rc=0.03;Lc=0.35e-3;
wnl=2*pi*50.5;rc=0.03;w0=2*pi*49.87;wo=w0;wn=wnl;

del1=0;del2=-0.011;

%Converter output voltage
Vod1=325.26;Voq1=0;Vod2=325.26;Voq2=0;

%Output voltage change
Ti1=inv([cos(del1) -sin(del1);sin(del1) cos(del1)]);
Ti2=inv([cos(del2) -sin(del2);sin(del2) cos(del2)]);


VoDQ1=Ti1*[[Vod1;Voq1]+[-Vod1*sin(del1)-Voq1*cos(del1);Vod1*cos(del1)-Voq1*sin(del1)]*[del1]];
VoDQ2=Ti2*[[Vod2;Voq2]+[-Vod2*sin(del2)-Voq2*cos(del2);Vod2*cos(del2)-Voq2*sin(del2)]*[del2]];

    
VoD1=VoDQ1(1,1);
VoQ1=VoDQ1(2,1);

VoD2=VoDQ2(1,1);
VoQ2=VoDQ2(2,1);


%COnverter output current
Iod1=39.21;Ioq1=-18.46;Iod2=70.66;Ioq2=-5.49;

%Output current change
IoDQ1=Ti1*[[Iod1;Ioq1]+[-Iod1*sin(del1)-Ioq1*cos(del1);-Iod1*cos(del1)-Ioq1*sin(del1)]*[del1]];
IoDQ2=Ti2*[[Iod2;Ioq2]+[-Iod2*sin(del2)-Ioq2*cos(del2);-Iod2*cos(del2)-Ioq2*sin(del2)]*[del2]];

IoD1=IoDQ1(1,1);
IoQ1=IoDQ1(2,1);

IoD2=IoDQ2(1,1);
IoQ2=IoDQ2(2,1);


%INVERTER OUTPUT Current
Ild1=39.00;Ilq1=-13.35;Ild2=70.44;Ilq2=-0.38;

%Inverter output reference change
IlDQ1=Ti1*[[Ild1;Ilq1]+[-Ild1*sin(del1)-Ilq1*cos(del1);Ild1*cos(del1)-Ilq1*sin(del1)]*[del1]];
IlDQ2=Ti2*[[Ild2;Ilq2]+[-Ild2*sin(del2)-Ilq2*cos(del2);Ild2*cos(del2)-Ilq2*sin(del2)]*[del2]];

IlD1=IlDQ1(1,1);
IlQ1=IlDQ1(2,1);

IlD2=IlDQ2(1,1);
IlQ2=IlDQ2(2,1);


%BUS Voltages
VbD1=322.14;Vbd2=322.58;VbQ1=-8.8;Vbq2=-12.65;

VbDQ2=Ti2*[[Vbd2;Vbq2]+[-Vbd2*sin(del2)-Vbq2*cos(del2);Vbd2*cos(del2)-Vbq2*sin(del2)]*[del2]];

VbD2=VbDQ2(1,1);
VbQ2=VbDQ2(2,1);


%LOAD Currents
Iloadd1=12.81; Iloadq1=-9.78;Iloadd2=20.70; Iloadq2=-13.65;Iloadd4=50.84; Iloadq4=-1.19;

%Inverter output reference change
Iload1DQ=Ti1*[[Iloadd1;Iloadq1]+[-Iloadd1*sin(del1)-Iloadq1*cos(del1);Iloadd1*cos(del1)-Iloadq1*sin(del1)]*[del1]];
Iload2DQ=Ti1*[[Iloadd2;Iloadq2]+[-Iloadd2*sin(del1)-Iloadq2*cos(del1);Iloadd2*cos(del1)-Iloadq2*sin(del1)]*[del1]];
Iload4DQ=Ti1*[[Iloadd4;Iloadq4]+[-Iloadd4*sin(del2)-Iloadq4*cos(del2);Iloadd4*cos(del2)-Iloadq4*sin(del2)]*[del2]];

IloadD1=Iload1DQ(1,1);
IloadQ1=Iload1DQ(2,1);

IloadD2=Iload2DQ(1,1);
IloadQ2=Iload2DQ(2,1);

IloadD4=Iload4DQ(1,1);
IloadQ4=Iload4DQ(2,1);

z1=zeros(3,2);
z2=zeros(2,2);
z3=zeros(7,2);
z4=zeros(10,1);
z5=zeros(1,10);
z6=zeros(2,10);
z7=zeros(2,13);
z8=zeros(13,2);
z9=zeros(1,13);
z10=zeros(1,10);
z11=zeros(13,13);

% GENERator 1 starts

AP1=[0 -mP1 0;0 -wc 0;0 0 -wc];
BP1=[0 0 0 0 0 0;0 0 wc*Iod1 wc*Ioq1 wc*Vod1 wc*Voq1;0 0 wc*Ioq1 -wc*Iod1 -wc*Voq1 wc*Vod1];
BPwcom1=[-1;0;0];
CPw1=[0 -mP1 0];
CPv1=[0 0 -nQ1;0 0 0];
BV11=[1 0;0 1];
BV21=[0 0 -1 0 0 0;0 0 0 -1 0 0];
CV1=[Kiv 0;0 Kiv];
DV11=[Kpv 0;0 Kpv];
DV21=[0 0 -Kpv -wn*Cf F 0;0 0 wn*Cf -Kpv 0 F];
BC11=[1 0;0 1];  
BC21=[-1 0 0 0 0 0;0 -1 0 0 0 0];
CC1=[Kii 0;0 Kii];
DC11=[Kpi 0;0 Kpi];
DC21=[-Kpi -wn*Lf 0 0 0 0;wn*Lf -Kpi 0 0 0 0];

ALCL1=[-rf/Lf w0 -1/Lf 0 0 0;-w0 -rf/Lf 0 -1/Lf 0 0;1/Cf 0 0 w0 -1/Cf 0;0 1/Cf -w0 0 0 -1/Cf;0 0 1/Lc 0 -rc/Lc w0;0 0 0 1/Lc -w0 -rc/Lc];
BLCL11=[1/Lf 0;0 1/Lf;0 0;0 0;0 0;0 0];
BLCL21=[0 0;0 0;0 0;0 0;-1/Lc 0;0 -1/Lc];
BLCL31=[Ilq1;-Ild1;Voq1;-Vod1;Ioq1;-Iod1];

TV1=[-VbD1*sin(del1)+VbQ1*cos(del1);-VbD1*cos(del1)-VbQ1*sin(del1)];
TS1=[cos(del1) -sin(del1);sin(del1) cos(del1)];
TC1=[-Iod1*sin(del1)-Ioq1*cos(del1);Iod1*cos(del1)-Ioq1*sin(del1)];

AGEN1=[AP1 z1 z1 BP1;BV11*CPv1 z2 z2 BV21;BC11*DV11*CPv1 BC11*CV1 z2 BC11*DV21+BC21;BLCL11*DC11*DV11*CPv1+BLCL21*[TV1 z2]+BLCL31*CPw1 BLCL11*DC11*CV1 BLCL11*CC1 ALCL1+BLCL11*(DC11*DV21+DC21)];
BGEN1=[z3;BLCL21*(inv(TS1))];
BwCOMM1=[BPwcom1;z4];
CGENw1=[CPw1 z5];
CGENc1=[TC1 z6 TS1];
CGEN1=[-Iod1*sin(del1)-Ioq1*cos(del1) z10 cos(del1) -sin(del1);Iod1*cos(del1)-Ioq1*sin(del1) z10 sin(del1) cos(del1)];

% GENERator 2

AP2=[0 -mP2 0;0 -wc 0;0 0 -wc];
BP2=[0 0 0 0 0 0;0 0 wc*Iod2 wc*Ioq2 wc*Vod2 wc*Voq2;0 0 wc*Ioq2 -wc*Iod2 -wc*Voq2 wc*Vod2];
BPwcom2=[-1;0;0];
CPw2=[0 -mP2 0];
CPv2=[0 0 -nQ2;0 0 0];
BV12=[1 0;0 1];
BV22=[0 0 -1 0 0 0;0 0 0 -1 0 0];
CV2=[Kiv 0;0 Kiv];
DV12=[Kpv 0;0 Kpv];
DV22=[0 0 -Kpv -wn*Cf F 0;0 0 wn*Cf -Kpv 0 F];
BC12=[1 0;0 1];  
BC22=[-1 0 0 0 0 0;0 -1 0 0 0 0];
CC2=[Kii 0;0 Kii];
DC12=[Kpi 0;0 Kpi];
DC22=[-Kpi -wn*Lf 0 0 0 0;wn*Lf -Kpi 0 0 0 0];

ALCL2=[-rf/Lf w0 -1/Lf 0 0 0;-w0 -rf/Lf 0 -1/Lf 0 0;1/Cf 0 0 w0 -1/Cf 0;0 1/Cf -w0 0 0 -1/Cf;0 0 1/Lc 0 -rc/Lc w0;0 0 0 1/Lc -w0 -rc/Lc];
BLCL12=[1/Lf 0;0 1/Lf;0 0;0 0;0 0;0 0];
BLCL22=[0 0;0 0;0 0;0 0;-1/Lc 0;0 -1/Lc];
BLCL32=[Ilq2;-Ild2;Voq2;-Vod2;Ioq2;-Iod2];

TV2=[-VbD2*sin(del2)+VbQ2*cos(del2);-VbD2*cos(del2)-VbQ2*sin(del2)];
TS2=[cos(del2) -sin(del2);sin(del2) cos(del2)];
TC2=[-Iod2*sin(del2)-Ioq2*cos(del2);Iod2*cos(del2)-Ioq2*sin(del2)];


AGEN2=[AP2 z1 z1 BP2;BV12*CPv2 z2 z2 BV22;BC12*DV12*CPv2 BC12*CV2 z2 BC12*DV22+BC22;BLCL12*DC12*DV12*CPv2+BLCL22*[TV2 z2]+BLCL32*CPw2 BLCL12*DC12*CV2 BLCL12*CC2 ALCL2+BLCL12*(DC12*DV22+DC22)];
BGEN2=[z3;BLCL22*(inv(TS2))];
BwCOMM2=[BPwcom2;z4];
CGENw2=z9;
CGENc2=[TC2 z6 TS2];
CGEN2=[-Iod2*sin(del2)-Ioq2*cos(del2) z10 cos(del2) -sin(del2);Iod2*cos(del2)-Ioq2*sin(del2) z10 sin(del2) cos(del2)];


AGEN=[AGEN1+BwCOMM1*CGENw1 z11;BwCOMM2*CGENw1 AGEN2];
BGEN=[BGEN1 z8;z8 BGEN2];
BwCOMM=[BwCOMM1;BwCOMM2];
CGENc=[CGENc1 z7;z7 CGENc2];
CGENw=[CGENw1 CGENw2];
CGEN=[CGENw;CGENc];
% generator Complete

% Transmission Lines
%LINE Currents
IlineD1=5.386;
IlineQ1=4.34;

%Line Parameters%
Rline1=0.23;
Lline1=0.35e-3;

% Line 1 

ANET1=[-Rline1/Lline1 wo;-wo -Rline1/Lline1];
B2NET1=[IlineQ1;-IlineD1];
B1NET1=[1/Lline1 0;0 1/Lline1];
CNET1=[1 0;0 1];

ANET=ANET1;
B1NET=[B1NET1 -1*B1NET1];
B2NET=B2NET1;
CNET=eye(2);


%  RL Load 
Rload1=16.21; Lload1=37.237e-3;
    
ALOAD1=[-Rload1/Lload1 w0;-w0 -Rload1/Lload1];
B1LOAD1=[1/Lload1 0 0 0;0 1/Lload1 0 0];
B2LOAD1=[IloadQ1;-IloadD1];
CLOAD1=[1 0;0 1];

% CPL Load 
Rload2=-13.4874; Lload2=-26.60e-3;
    
ALOAD2=[-Rload2/Lload2 w0;-w0 -Rload2/Lload2];
B1LOAD2=[1/Lload2 0 0 0;0 1/Lload2 0 0];
B2LOAD2=[IloadQ2;-IloadD2];
CLOAD2=[1 0;0 1];

% RIAL LOAD

KpvAL=0.05; KivAL=0.5; 
KpiAL=0.2023; KiiAL=150;

LfAL=2.3e-3; Cdc=0.1e-6;rdc=1*10000*1000; LcAL=0.93e-3;rcAL=0.03;
CfAL=8.8e-6;rfAL=0.1;
IlDAL=25.52; IlQAL=0.01;
VoDAL=322.58; VoQAL=-4.98;
ViDAL=325.26; ViQAL=-12.65;
IoDAL=25.12; IoQAL=0.909;
Vdc=700;
delAL=-0.322;

TiAL=inv([cos(delAL) -sin(delAL);sin(delAL) cos(delAL)]);

% Vo REF CHANGE RIAL
VodqAL=TiAL*[VoDAL;VoQAL];
    
VodAL=VodqAL(1,1);
VoqAL=VodqAL(2,1);

% Il REF CHANGE RIAL
IldqAL=TiAL*[IlDAL;IlQAL];
    
IldAL=IldqAL(1,1);
IlqAL=IldqAL(2,1);

% Io REF CHANGE RIAL
IodqAL=TiAL*[IoDAL;IoQAL];
    
IodAL=IodqAL(1,1);
IoqAL=IodqAL(2,1);

% Vo REF CHANGE RIAL
VidqAL=TiAL*[ViDAL;ViQAL];
    
VidAL=VidqAL(1,1);
ViqAL=VidqAL(2,1);


ARIAL11=zeros(1,9);
ARIAL14=-1;
ARIAL21=[1 0;0 1]*[KivAL;0];
ARIAL22=zeros(2,2);
ARIAL23=[-1 0 0 0 0 0;0 -1 0 0 0 0];
ARIAL24=[-KpvAL;0];
ARIAL31=[KivAL*KpiAL;zeros(5,1)];
ARIAL32=[KiiAL*[1 0;0 1]/LfAL zeros(2,4)]';
ARIAL33=[(-rfAL-KpiAL)/LfAL w0-wn 1/LfAL 0 0 0;wn-w0 (-rfAL-KpiAL)/LfAL 0 1/LfAL 0 0;-1/CfAL 0 0 w0 1/CfAL 0;0 -1/CfAL -w0 0 0 1/CfAL;0 0 -1/LcAL 0 -rcAL/LcAL w0;0 0 0 -1/LcAL -w0 -rcAL/LcAL];
ARIAL34=[KpvAL*KpiAL/LfAL zeros(1,5)]';
ARIAL41=[-IldAL*KivAL*KpvAL/(Vdc*Cdc)];
ARIAL42=[-IldAL*KiiAL/(Vdc*Cdc) -IlqAL*KiiAL/(Vdc*Cdc)];
ARIAL43=[(VidAL+IldAL*KpiAL-IlqAL*LfAL*wn)/(Vdc*Cdc) (ViqAL+IlqAL*KpiAL+IldAL*LfAL*wn)/(Vdc*Cdc) zeros(1,4)];
ARIAL44=[(IldAL*KpiAL*KpvAL/(Vdc*Cdc))-((IldAL*VidAL+IlqAL*ViqAL)/(Vdc*Vdc*Cdc))-1/(rdc*Cdc)];


AAL=[ARIAL11 ARIAL14;ARIAL21 ARIAL22 ARIAL23 ARIAL24;ARIAL31 ARIAL32 ARIAL33 ARIAL34;ARIAL41 ARIAL42 ARIAL43 ARIAL44];

BALw=[zeros(1,3) IlqAL -IldAL VoqAL -VodAL IoqAL -IodAL 0]';
BALv=[zeros(1,7) (-cos(delAL))/LcAL (sin(delAL))/LcAL 0;zeros(1,7) (-sin(delAL))/LcAL (-cos(delAL))/LcAL 0;zeros(2,10)]';
BALu=[1 KpvAL 0 KpiAL*KpvAL/LfAL 0 zeros(1,4) (-IldAL*KpiAL*KpvAL)/(Vdc*Cdc);0 0 1 0 KpiAL/LfAL zeros(1,4) (-IlqAL*KpiAL)/(Vdc*Cdc)]';

CAL=[zeros(1,7) cos(delAL) -sin(delAL) 0;zeros(1,7) sin(delAL) cos(delAL) 0];

% R Load
Rload4=6.348; Lload4=1.0e-6;
    
ALOAD4=[-Rload4/Lload4 w0;-w0 -Rload4/Lload4];
B1LOAD4=[0 0 1/Lload4 0 ;0 0 0 1/Lload4];
B2LOAD4=[IloadQ4;-IloadD4];
CLOAD4=[1 0;0 1];

APL=[ALOAD1 z2 z2;z2 ALOAD2 z2;z2 z2 ALOAD4];
B1PL=[B1LOAD1;B1LOAD2;B1LOAD4];
B2PL=[B2LOAD1;B2LOAD2;B2LOAD4];
CPL=[CLOAD1 z2 z2;z2 CLOAD2 z2;z2 z2 CLOAD4];

% Load Ends

% Mapping matrix
I1=eye(4);
RN=1000*I1;
MGEN=I1;
MNET=[-1 0;0 -1;1 0;0 1];
MPL=[-1 0 -1 0 0 0;0 -1 0 -1 0 0;0 0 0 0 -1 0;0 0 0 0 0 -1];
MAL=[0 0;0 0;-1 0;0 -1];
MDIST=[-1 0;0 -1;0 0;0 0];
CDIST=[1 0;0 1];

AMG=[AGEN+BGEN*RN*MGEN*CGENc+BwCOMM*CGENw BGEN*RN*MNET*CNET BGEN*RN*MPL*CPL BGEN*RN*MAL*CAL;B1NET*RN*MGEN*CGENc+B2NET*CGENw ANET+B1NET*RN*MNET*CNET B1NET*RN*MPL*CPL B1NET*RN*MAL*CAL;B1PL*RN*MGEN*CGENc+B2PL*CGENw B1PL*RN*MNET*CNET APL+B1PL*RN*MPL*CPL B1PL*RN*MAL*CAL;BALv*RN*MGEN*CGENc+BALw*CGENw BALv*RN*MNET*CNET BALv*RN*MPL*CPL AAL+BALv*RN*MAL*CAL];
BMG=[BGEN*RN*MDIST zeros(26,2);B1NET*RN*MDIST zeros(2,2);B1PL*RN*MDIST zeros(6,2);BALv*MDIST*CDIST BALu];
CMG=[CGEN zeros(5,2) zeros(5,6) zeros(5,10);zeros(6,26) zeros(6,2) CPL zeros(6,10);zeros(2,26) zeros(2,2) zeros(2,6) CAL;MGEN*CGENc MNET*CNET MPL*CPL MAL*CAL];
DMG=[zeros(13,4);MDIST*CDIST  zeros(4,2)];

size(MDIST*CDIST);
X=eig(AMG);
plot(X,'*');
grid on;


[V,D] = eig(AMG);
[d,ind] = sort(diag(D));
x=real(d);
y=imag(d);


xlswrite('PoleZeroAnalysis.xlsx',[x(:),y(:),ind]);  % Making columns of real and imaginary value of states and state number










