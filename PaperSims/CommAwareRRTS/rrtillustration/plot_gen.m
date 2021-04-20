mrksz = 20;
p0 = [1,1];
p1=[0.9, 2];
p2=[2, 2.5];

p3 = [2.2, 1.4];
p4 = p3+[0.4, 0.5];
p5 = p4+[1,0];
p6 = p3+[1,0.1];
plot([1,0.9], [1, 2],'.-k', 'markersize',mrksz);
hold on
plot([0.9, 1.8], [2,2.5],'.-k', 'markersize',mrksz);

plot([1, 2.2],[1,1.4],'.-k','markersize',mrksz);
plot([2.2, p4(1)],[1.4, p4(2)],'.-k', 'markersize',mrksz);
plot([p4(1), p5(1)],[p4(2), p5(2)],'.-k', 'markersize',mrksz);
plot([p3(1), p6(1)],[p3(2), p6(2)],'.-k', 'markersize',mrksz);

%Add for frame 2
%%
p_new = [3,1.7];
plot(p_new(1), p_new(2),'.b','markersize',mrksz)
%%
% Frame 3
plot([p3(1), p_new(1)], [p3(2), p_new(2)], 'b');
%%
% Frame 4
theta = 0:pi/50:2*pi;
d = 0.9;
x = d*cos(theta)+p_new(1);
y = d*sin(theta)+p_new(2);
plot(x,y,'g');

