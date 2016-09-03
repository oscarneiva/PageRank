%PageRank Simulation (Matlab version)
%Laboratório Nacional de Computação Científica - LNCC
%Oscar Neiva Eulálio Neto - Orientador: Marcos Garcia Todorov

clear
clc
close all

%Start get the running time in seconds
tic();

%Number of states
N = 4;

%Number of distributed matrix
NPn = N;

%Time horizon
T = 50;

%Time horizon monte carlo 
itmax = 50;

%Parameter m of Mean-Square Error
m = 0.15; %m = 0.15

%Random vector 
%x = rand(1,N);
%x = diag(x*ones(N,1))\x;

%Static vectors
x = [0 0 0 1]; %Simple PageRank

%Other vectors
y = x;
xd = x;
yd = x;
z = x;
zr = x;
w = zeros(T+1,N);

%Random Stochastic Matrix(line)
%P = rand(N);
%P = diag(P*ones(N,1))\P;

%Stochastic Matrix
P = [0 0 0 1; 0.5 0 0.5 0; 0 0 0 1; 0.5 0.5 0 0];

%Vector with ones
v1 = ones(1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mean-Square Error Convergence
M = 2*m/(N-m*(N-2));

%Create a tridimentional matrix with zeros
PP = zeros(size(P,1),size(P,2),N);

%Insert values in the tridimensional matrix
for i=1:N
	%Insert ones in the diagonals
	PP(i,i,:) = 1;
	
	%Copy lines to the distributed matrix
    PP(i,:,i) = P(i,:);
end

for it=1:itmax
	for k=1:T
		%Generate one integer number to represent the random distributed matrix
		Pn = ceil(rand*NPn);

		%Power Metod
		x(k+1,:) = x(k,:)*P;
		
		%PageRank with teleportation model (not distributed)
		y(k+1,:) = (1-M)*(x(k,:)*P) + (M/N)*v1;

		%PageRank with distributed link matrices
		xd(k+1,:) = xd(k,:) * PP(:,:,Pn);

		%PageRank with distributed teleportation model
		yd(k+1,:) = (1-M)*(xd(k,:)*PP(:,:,Pn)) + (M/NPn)*v1;

		%PageRank with time average (not recursive)
		z(k,:) = sum(yd)/(T+1);

		%PageRank with time average (recursive)
		zr(k+1,:) = (((k+1)/(k+2)) * zr(k,:)) + ((1/(k+2))*yd(k+1,:));
	end
	w = w + zr./itmax;
end


time = toc ();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot results
clf

figure(1);
%subplot(211);
plot(x,'Linewidth',2.0);
set(gca,'fontsize',20);
set(gca,'XLim',[1 T]);
set(gca,'YLim',[0 2/N]);
title('Power Method');
ylabel('Pages values');
xlabel('Time');
print('-depsc2','powermetod');
%!ps2pdf -dEPSCrop powermetod.eps powermetod.pdf

figure(2);
%subplot(212);
plot(y, 'Linewidth', 2.0);
set(gca,'fontsize',20);
set(gca,'XLim',[1 T]);
set(gca,'YLim',[0 2/N]);
title('Teleportation Model');
ylabel('Pages values');
xlabel('Time');
print('-depsc2','teleportation');
%!ps2pdf -dEPSCrop teleportation.eps teleportation.pdf

%figure(3);
%subplot(211);
%plot(xd, 'Linewidth', 1.0);
%set(gca,'fontsize',20);
%set(gca,'XLim',[1 T]);
%set(gca,'YLim',[0 2/N]);
%title('PageRank With Distributed Link Matrices');
%ylabel('Pages values');
%xlabel('Time');
%print('-depsc2','pagedistributed');
%!ps2pdf -dEPSCrop pagedistributed.eps pagedistributed.pdf

figure(4);
%subplot(212);
plot(yd, 'Linewidth', 2.0);
set(gca,'fontsize',20);
set(gca,'XLim',[1 T]);
set(gca,'YLim',[0 2/N]);
title('Teleportation Model With Distributed Link Matrices');
ylabel('Pages values');
xlabel('Time');
print('-depsc2','teledistributed');
%!ps2pdf -dEPSCrop teledistributed.eps teledistributed.pdf

%figure(5);
%subplot(7,2,9);
%plot(z);
%set(gca,'fontsize',20);
%set(gca,'XLim',[1 T]);
%set(gca,'YLim',[0 2/N]);
%title('Time Average (Not Recursive)');
%ylabel('Pages values');
%xlabel('Time');
%print('-depsc2','polyak','-F:30');
%!ps2pdf -dEPSCrop polyak.eps polyak.pdf

figure(6);
%subplot(211);
plot(zr, 'Linewidth', 2.0);
set(gca,'fontsize',20);
set(gca,'XLim',[1 T]);
set(gca,'YLim',[0 2/N]);
title('Time Average (Recursive)');
ylabel('Pages values');
xlabel('Time');
print('-depsc2','timerecursive');
%!ps2pdf -dEPSCrop timerecursive.eps timerecursive.pdf

figure(7);
%subplot(212);
plot(w, 'Linewidth', 2.0);
set(gca,'fontsize',20);
set(gca,'XLim',[1 itmax]);
set(gca,'YLim',[0 2/N]);
title('Monte Carlo Simulation');
ylabel('Pages values');
xlabel('Time');
print('-depsc2','montecarlo');
%!ps2pdf -dEPSCrop montecarlo.eps montecarlo.pdf
