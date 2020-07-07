function dx = dyn_ballisticM_thrust(x,b,thr, mass) 
    global t bal_inform accBM BC Thr totalMass
    
    %---load earth information---------------------------------------------
    earth_constants      ;
    
    %dM/dt of ballistic missile
    dM.Stage1 = bal_inform.launchinform.FuWeiStg1 / bal_inform.launchinform.BurTimStg1    ;
    dM.Stage2 = bal_inform.launchinform.FuWeiStg2 / bal_inform.launchinform.BurTimStg2    ;

    %canister weight of bal.listic missile
    bal_inform.canWeiStg1 = bal_inform.launchinform.GrWeiStg1 - bal_inform.launchinform.FuWeiStg1    ;
    bal_inform.canWeiStg2 = bal_inform.launchinform.GrWeiStg2 - bal_inform.launchinform.FuWeiStg2    ;

    %next stage time
    Stage_Time.one   = bal_inform.launchinform.BurTimStg1     ;
    Stage_Time.two = 0                                        ; %free flight
    %----------------------------------------------------------------------
    
    
    %---update ballistic information (total mass & thrust)-----------------
    if t <= Stage_Time.one

            ISPT = bal_inform.launchinform.ISPstg1   ;
            dMdt = dM.Stage1            ;            
    else
                ISPT = bal_inform.launchinform.ISPstg2                  ;
                dMdt = dM.Stage2                                        ;
    end
    
    totMass = totalMass;
    magThrBM = Thr;    

    magPosBM = sqrt(x(1)^2 + x(2)^2 + x(3)^2)           ; % magnitude of position vector of ballistic missile
    unPosBM = x(1:3)./ magPosBM                         ; % unit vector of bal.listic missile position vector
    gBM = (constant.Gc * constant.Me) / (magPosBM ^ 2)  ; % gravitational acceleration of BM
    magVelBM = sqrt(x(4)^2 + x(5)^2 + x(6)^2)           ; % magnitude of velocity vector of BM
    unBMvel = x(4:6)./magVelBM                          ; % unit vector of vel vec of BM
    magWeiBM = totMass * gBM                            ; % magnitude of weight vector of ballistic missile
    unMagWeiBM = -unPosBM                               ; % unit vec of weight vector of ballistic missile
    weiVec = unMagWeiBM * magWeiBM                      ; % weight vector of ballistic missile
%     magThrBM = dMdt * gBM * ISPT                        ; % magnitude of thrust vector of ballistic missile  
    unThrBM = unBMvel                                   ; % unit vector of thrust vec of ballistic missile   
    AltitudeBM = norm(x(1:3))-constant.Re               ; % altitude of ballistic missile
    %----------------------------------------------------------------------
    
    
    %---drag define--------------------------------------------------------
    if t<= bal_inform.launchinform.BurTimStg1
        tableCd = bal_inform.Cd.Cd1 ;
        d = bal_inform.Cd.diameter(1); 
    else
        tableCd = bal_inform.Cd.Cd2 ;
        d = bal_inform.Cd.diameter(2);
    end
    %----------------------------------------------------------------------
    
    
    %---Earth Atmosphere Model---------------------------------------------
    if AltitudeBM >= 25000 
        T = -131.21 + 0.00299*AltitudeBM        ;
        P = 2.488*((T + 273.1)/216.6)^-11.388   ;
        elseif (AltitudeBM >= 11000) && (AltitudeBM < 25000) 
            T = -56.46                              ;
            P = 22.65*exp(1.73-0.000157*AltitudeBM) ;
            else
                T = 15.04 - 0.00649*AltitudeBM          ;
                P = 101.29*((T + 273.1)/288.08)^5.256   ;
    end
    rho = P/(0.2869*(T+273.1)) ;
    Tk = 273.1+T  ;
    speedS = sqrt(286*Tk*(1 + (1.4-1))/(1 + (1.4-1)*((3056/Tk)^2*(exp(3056/Tk)/(exp(3056/Tk)-1)^2)))) ; % speed of sound
    machN = magVelBM/speedS   ;  % mach number
    tempCd = [];
    if AltitudeBM <= 85000    
        for i=1:1:size(bal_inform.Cd.mach,2)
            tempCd = [tempCd interp1(bal_inform.Cd.alt,tableCd(:,i),AltitudeBM)];
        end
        if machN > bal_inform.Cd.mach(end)
            machN = bal_inform.Cd.mach(end);
        end
        Cd = interp1(bal_inform.Cd.mach,tempCd,machN) ;
        if BC < 0
            BMcoeff = totMass/(Cd*(pi*d^2/4)); % ballistic coefficient  
            accDragBM = rho*magVelBM^2/(2*BMcoeff)*unThrBM;    
        elseif BC >= 0
            BMcoeff = BC;
            accDragBM = (1/2)*rho*magVelBM^2/(BMcoeff)*unThrBM;
        end
    else
        accDragBM = 0;
    end
    thrBM  = magThrBM * unThrBM; % thrust vector of ballistic missil  
    coriBM = -cross(constant.omega,x(4:6))...     ; % coriolis acceleration    
             -cross(constant.omega,cross(constant.omega,x(1:3))) ;
    %---calculate acceleration of ballistic missile------------------------
    accBM = (weiVec + thrBM)/totMass +coriBM  -  accDragBM  ; % acceleration    
    dx(1) = x(4)    ;
    dx(2) = x(5)    ;
    dx(3) = x(6)    ;  
    dx(4:6) = accBM ;
    dx = [dx(1) dx(2) dx(3) dx(4) dx(5) dx(6)]';
return
