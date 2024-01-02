clear
clc
%% Ex 1
n = [4,5,4,5];
m = [0,5,2,4];
theta = (0:180);
for i = 1:length(theta)
    for j = 1:4
    Pnm(j,i) = legendre(n(j),m(j),theta(i));
    end
end

figure (1)
hold on
plot(Pnm(1,:),theta,'-r')
plot(Pnm(2,:),theta,'-g')
plot(Pnm(3,:),theta,'-b')
plot(Pnm(4,:),theta,'-.m','LineWidth',2)
set(gca,'Ydir','reverse')
xlabel('value of LEGENDRE function')
ylabel('colatitude [deg]')
legend('n4m0','n5m5','n4m2','n5m1')
hold off

