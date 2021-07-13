clc
close all
clear all

Ts = 0.1 ;	% [s]
Tend  = 10 ;	% [s]
N  = Tend / Ts ;

sys.A = [ -0.5572  -0.7814 ;
	   0.7814   0      ];
sys.B = [ 1  -1 ;
	  0   2 ];
sys.C = [ 1.9691  6.4493 ];
sys.D = [ 0  0 ];

u0_h = zeros(2,N);
u0_h(:,1) = [ 1 ; 1 ];
%u0_h = randn(2,N);
x0   = [ 0; 0 ];
y0_h = drivePlant(sys, x0, u0_h);

s = 20 ;
m = 40 ;

H  = [];
Hd = [];

for i = 1:s
	H  = [ H  ; y0_h(i:i+m-s-1) ];
	Hd = [ Hd ; y0_h(i+1:i+m-s) ];
end

% PDQ' , USV', 
[  U,  S,  V ] = svd(H);
[ Ud, Sd, Vd ] = svd(Hd);

r  = rank(H);
rd = rank(Hd);

Ur = U(:,1:r);
Vr = V(:,1:r);
Sr = S(1:r,1:r);


%Ar = pinv(S).^0.5 * U' * Hd * V * pinv(S).^0.5
mod.A = pinv(Sr).^0.5 * Ur' * Hd * Vr * pinv(Sr).^0.5
mod.B = Sr.^0.5 * Vr' * [ eye(2) ; zeros(m-s-2,2) ]
mod.C = [ eye(1)  zeros(1,m-s-1) ] * Ur * Sr.^0.5
mod.D = mod.C * mod.B

y_h = drivePlant(mod, zeros(4,1), u0_h);

figure(1)
hold on
plot(y_h);
plot(y0_h,':');
hold off

function y_h = drivePlant(sys, x0, u_h)
	N = length(u_h) ;
	x = x0 ;
	x_h = zeros(length(x),N) ;
	for i=1:N-1
		x1 = sys.A * x + sys.B * u_h(:,i) ;
		y  = sys.C * x + sys.D * u_h(:,i) ;

		y_h(:,i+1) = y;
		x = x1 ;
	end
end
