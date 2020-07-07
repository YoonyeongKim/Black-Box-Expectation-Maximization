function dx = dyn_ballisticM_HV61_thrust_baseline(x,time,deltaTime,ballistic_coefficient, thrust_magnitude, mass)
    global idxMissile curMissile
    global t data_name
    global bal_inform accBM BC Thr totalMass
    persistent BM 
    
    if isempty(BM) 
        BM = load(data_name); % not efficient
    end
    
    BC = ballistic_coefficient;
    Thr = thrust_magnitude;
    totalMass = mass;
    t = time;
    earth_constants ; % load earth information

    % site1 to seoul azimuth angle: 60 / site4 to seoul azimuth angle: 128
    RadarPos = [37.53 126.96 94 10] ; %latitude(degree), longitude(degree), azimuth(degree), altitude(m)
    
    [th,ph] = geo2sph('N',RadarPos(1),'E',RadarPos(2))      ;
    [Rx,Ry,Rz] = sph2car(th,ph,constant.Re)                 ;
    [RadarPosX,RadarPosY,RadarPosZ] = ecef2enu(RadarPos(1)*pi/180,RadarPos(2)*pi/180,RadarPos(3)*pi/180,Rx,Ry,Rz) ;  

    x(1:3,1) = x(1:3,1) + [RadarPosX;RadarPosY;RadarPosZ + RadarPos(4)];
    x(4:9,1) = x(4:9,1);

    x(1:3,1) = enu2ecef(RadarPos(1)/180*pi,RadarPos(2)/180*pi,RadarPos(3)/180*pi,x(1,1),x(2,1),x(3,1));
    x(4:6,1) = enu2ecef(RadarPos(1)/180*pi,RadarPos(2)/180*pi,RadarPos(3)/180*pi,x(4,1),x(5,1),x(6,1));
    x(7:9,1) = enu2ecef(RadarPos(1)/180*pi,RadarPos(2)/180*pi,RadarPos(3)/180*pi,x(7,1),x(8,1),x(9,1));   
    
    selected_target_inform = [BM.ballistic_target(idxMissile(curMissile)).BMtype; BM.ballistic_target(idxMissile(curMissile)).launchSiteid; BM.ballistic_target(idxMissile(curMissile)).destinationid; BM.ballistic_target(idxMissile(curMissile)).time(1)];
    target = BMinform(1,selected_target_inform) ;
    bal_inform =  target(1) ;   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp = RK4_thrust(@dyn_ballisticM_thrust, x(1:6,1), 0, 0, 0, deltaTime);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    temp(7:9,1)= accBM;
   
    [temp(1,1), temp(2,1), temp(3,1)] = ecef2enu(RadarPos(1)*pi/180,RadarPos(2)*pi/180,RadarPos(3)*pi/180,temp(1,1), temp(2,1), temp(3,1)) ;
    [temp(4,1), temp(5,1), temp(6,1)] = ecef2enu(RadarPos(1)*pi/180,RadarPos(2)*pi/180,RadarPos(3)*pi/180,temp(4,1), temp(5,1), temp(6,1)) ;
    [temp(7,1), temp(8,1), temp(9,1)] = ecef2enu(RadarPos(1)*pi/180,RadarPos(2)*pi/180,RadarPos(3)*pi/180,temp(7,1), temp(8,1), temp(9,1)) ;

    temp(1,:) = temp(1,:) - RadarPosX ;
    temp(2,:) = temp(2,:) - RadarPosY ;
    temp(3,:) = temp(3,:) - RadarPosZ - RadarPos(4) ; % consider altitude of radar    
     
    dx = temp;
end