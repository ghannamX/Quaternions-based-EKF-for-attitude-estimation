%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D orientation estimation template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONFIGURATION and INITIALISATION
dt=0.01;
COM='COM5';

useSimulation = true;
useSimulation_visuals = useSimulation;
imu411dosimulate = false;
%% CALIBRATION
% To re-run calibration, just clear MAG0 or run the calibration script
% separately:
if ~exist('MAG0', 'var')
    runcalibration;
end
%% affichage du parallélépipède
launch3D;
% Fairly Good Tuning : 
% qa = 0.001    ; ra = 1000 accvar ; rm = 1e+4 magvar  ; P = 1000*eye(4);
% qa =  0.002  ; ra = 1000 accvar ; rm = 1e+4 magvar  ; P = 1000*eye(4);
%% estimation attitude
qa = 0.001         ;     %Bruit d'état
ra = 1000 * ACCVAR ;     %Bruit de mesure accelero
rm = 1e+4 * MAGVAR;      %Bruit de mesure magnetomter
Q = diag([qa qa qa qa]); %
R = diag([ra' rm']);
X=[1 0 0 0]';            %Etat initial : q
P = 1000*eye(4);                                                    

tp=0;
ii=0;
obs = [];
xtrue = [];
xhat = [];  
if useSimulation
    imu411('open',COM, 'simimu_3Dmaneuver');
else
    imu411('open',COM);
end

while(1)
    ii=ii+1;
    taxis = 1:ii;
    % read sensor
    [d, xt] = imu411('read'); %cette lecture est bloquante et cadence la boucle a 0.01s
    obs(:, end+1) = d;
    xtrue(:, end+1) = xt ;   %Rappel
    xt;
    % d(1)    = Time                (s)
    % d(2:4)  = Gyroscope     X,Y,Z (°/s)
    % d(5:7)  = Accelerometer X,Y,Z (g)
    % d(8:10) = Magnetometer  X,Y,Z (Gauss)
    t=d(1);
    % Predict
    X = X;
    F = eye(4)  ;
    P = F*P*F'+Q;
    % Update
    if ~isnan(t)
        Y =[d(5)    %accX
            d(6)    %accY
            d(7)    %accZ
            d(8)    %magX
            d(9)    %magY
            d(10)]; %magZ


        % size(Y) = 6 * 1
        M=quat2M(X(1:4))';
        Yhat = [M*ACC0
            M * MAG0];
        q0=X(1);
        q1=X(2);
        q2=X(3);
        q3=X(4);
        J1=2*[q0 q1 -q2 -q3
            -q3 q2 q1 -q0
            q2 q3 q0 q1];
        J2=2*[q3 q2 q1 q0
            q0 -q1 q2 -q3
            -q1 -q0 q3 q2];
        J3=[-q2 q3 -q0 q1
            q1 q0 q3 q2
            q0 -q1 -q2 q3];
        Hacc=J1*ACC0(1)+J2*ACC0(2)+J3*ACC0(3);
        Hmag =  J1*MAG0(1)+J2*MAG0(2)+J3*MAG0(3); % M * MAG0% The magnetometer measures                                                    
        H=[Hacc     % , in the moving frame of the senso 4 * 3
           Hmag]; % size = 6 * 4
        G = H*P*H'+R;
        K = P*H'*inv(G);

        X = X+K*(Y-Yhat);
        X(1:4) = X(1:4)/norm(X(1:4));
        P = P - P*K*H;
        xhat(:, end+1) = X;        
    end

    DCM_k = quat2dcm(X(1:4)'/norm(X(1:4)))';
    [yawi(ii), pitchi(ii), rolli(ii)] = quat2angle(X(1:4)'/norm(X(1:4)) ); % return Yaw and


end



