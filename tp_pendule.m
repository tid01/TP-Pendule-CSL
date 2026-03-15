%Definition des variables du systeme

mr=0.095;
Lr=0.085;
Jr=5.72e-5;
mp=0.024;
Lp=0.129;
Jp=3.33e-5;
kt=0.042;
km=0.042;
Rm=8.4;
g=9.81;

Jt=Jp*Jr+Jr*mp*(Lp^2)/4+mp*Jp*Lr^2;


%Definition des matrices du systeme a l'aide des equations donnees

M1=[0 ((mp^2)*(Lp^2)*Lr*g)/(4*Jt); 0 (mp*Lp*g*(mp*(Lr^2)+Jr))/(2*Jt)];
M2=[-(kt*km*(Jp+mp*(Lp^2)/4))/(Jt*Rm) 0; -(kt*km*mp*Lp*Lr)/(2*Rm*Jt) 0];
M3=[(kt*(Jp+mp*(Lp^2)/4))/(Rm*Jt);(kt*mp*Lp*Lr)/(2*Rm*Jt)];

A=[0 0 1 0;0 0 0 1; M1 M2];

B=[0;0;M3];
C = [1 0 0 0];
D = 0;

%Definition du systeme
sys=ss(A,B,C,D);

%Verifications du systeme
SpA = eig(A); %verif stabilite

Qc=ctrb(A,B); %verif commandabilite



%Changement de variable x=M*z pour FCC
p=charpoly(A);
[a3, a2,a1,a0]=deal(p(2),p(3),p(4),p(5));
m4=B;
m3=(A+a3*eye(4))*B;
m2=(A^2+a3*A+a2*eye(4))*B;
m1=(A^3+a3*A^2+a2*A+a1*eye(4))*B;

%M=[m1 m2 m3 m4];   %les calculs renvoyes par matlab mettent des tres petit
%chiffres au lieu de 0, donc redefinition manuelle des matrices
M=[-5670.53242754656	0	49.7178109268685	0;
    0	0	49.1330536535634	0;
    0	-5670.53242754656	0	49.7178109268685;
    0	0	0	49.1330536535634];
Minv=inv(M);



%Az=Minv*A*M;
Az=[0	1	0	0;
0	0	1	0;
0	0	0	1;
0	238.162361956956	261.524955558447	-2.08814805892848];
%Bz=Minv*B;
Bz=[0;0;0;1];
Cz=C*M;
Dz=D;


%Definition du systeme sous forme compagne de commande
sys_z=ss(Az,Bz,Cz,Dz);



%Placement des poles pour un retour d'etat
Kz=place(Az,Bz,[-12,-13,-14,-15]);
K=Kz*Minv;