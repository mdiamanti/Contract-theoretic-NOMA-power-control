function [] = scatter_plot(femto,x_femto,y_femto,z_femto,users,x_users,y_users,z_users,chosen_femto,sla_ite)

color = [0 0.4470 0.7410 ; 0.8500 0.3250 0.0980 ; 0.9290 0.6940 0.1250 ; 0.4940 0.1840 0.5560 ; 0.4660 0.6740 0.1880 ; 0.6350 0.0780 0.1840];

figure();
for i = 1:femto
    scatter3(x_femto(i),y_femto(i),z_femto(i),80,'filled','d','MarkerEdgeColor',color(i,:),...
              'MarkerFaceColor',color(i,:),'LineWidth',1.5);
    hold on;
end

hold on;

for i = 1:users
    if chosen_femto(i) == 1
       scatter3(x_users(i),y_users(i),z_users(i),40,'filled','p','MarkerEdgeColor',color(1,:),'MarkerFaceColor',color(1,:));
    elseif chosen_femto(i) == 2
       scatter3(x_users(i),y_users(i),z_users(i),40,'filled','p','MarkerEdgeColor',color(2,:),'MarkerFaceColor',color(2,:));
    elseif chosen_femto(i) == 3
       scatter3(x_users(i),y_users(i),z_users(i),40,'filled','p','MarkerEdgeColor',color(3,:),'MarkerFaceColor',color(3,:)); 
    elseif chosen_femto(i) == 4
       scatter3(x_users(i),y_users(i),z_users(i),40,'filled','p','MarkerEdgeColor',color(4,:),'MarkerFaceColor',color(4,:));
    elseif chosen_femto(i) == 5
       scatter3(x_users(i),y_users(i),z_users(i),40,'filled','p','MarkerEdgeColor',color(5,:),'MarkerFaceColor',color(5,:));
    end
    hold on;
end

hold on;

xlabel('x axis [m]');
ylabel('y axis [m]');
zlabel('z axis [m]');
set(gca,'FontSize',20);

if sla_ite ~= -1
    title(['User to BS association - SLA iteration ',num2str(sla_ite)])
else
    title('User to BS association')
end

end

