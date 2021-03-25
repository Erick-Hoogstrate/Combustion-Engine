function v = engine_kinematics(bore,stroke,r,rod,start_crank,end_crank)
    a=stroke/2;
    R=rod/a;
    theta=linspace(start_crank,end_crank,abs(end_crank-start_crank));
    v_s = (pi/4)*bore^2*stroke;
    v_c=v_s/(r-1);
    t1=0.5*(r-1);
    t2=R+1-cosd(theta);
    t3=(R^2-sind(theta).^2).^0.5;
    v=(1+t1*(t2-t3)).*v_c;
end