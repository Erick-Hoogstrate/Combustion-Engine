figure()
hold on
plot(Volume, average_pressure)
plot(v0,p0/bara,'color','r')
plot(v1,p1/bara,'color','r')
plot(v2,p2/bara,'color','r')
plot(v3,p3/bara,'color','r')
plot(v4,p4/bara,'color','r')
set(gca,'FontSize',30)
P0=plot([v0 v1],[p0 p1]/bara,'color','red');%isobaric expansion
P1=plot(v_compression,p_compression/bara,'color','red');%isentropic compression
P2=plot([v2 v3],[p2 p3]/bara,'color','red');%isochoric heat addition
P3=plot(v_expansion,p_expansion/bara,'color','red');%isentropic expansion
P4=plot([v4 v1],[p4,p1]/bara,'color','red');%isochoric heat rejection
P5=plot([v1 v0],[p1 p0]/bara,'color','red');%isobaric compresion
xlabel('Volume [m^3]')
ylabel('Pressure [Bar]')
legend('Average','Otto')
title(['Combined pV Diagram E15 Full Load'])