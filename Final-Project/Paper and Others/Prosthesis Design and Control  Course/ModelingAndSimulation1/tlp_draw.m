function tlp_draw(x,t)
    global tlp
	
	L1 = tlp.L1;
	L2 = tlp.L2;
	
	q1 = x(1);
	q2 = x(2);
	plot([0 L1*cos(q1) L1*cos(q1)+L2*cos(q1+q2)], [0 L1*sin(q1) L1*sin(q1)+L2*sin(q1+q2)], 'k-o');
    hold on
    plot([0 0.15 -0.07 0],[0 -0.07 -0.07 0],'k');       % draw a foot (for visuals only)
    hold off
	axis('equal')
    text(0.6,0.9,['t=' num2str(t)],'FontSize',12);
	axis([-1 1 -1 1]);
    drawnow;
end