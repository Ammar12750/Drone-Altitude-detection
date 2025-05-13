

% MTRN4010.2025,
% Example that shows how to read the synthetic data for testing your
% solution to Project2/partA.1

function main()

r=GenerateSyntheticDataProject2_A1();        % get data

B = EstimateBiases(r);              % run estimation.

% print your result (estimated biases)
disp('---------------------------------------')
fprintf('Est. biases =[%.2f][%.2f][%.2f] deg/sec. (Estimated by optimizer)\n',B*180/pi);
fprintf('Real biases =[%.2f][%.2f][%.2f] deg/sec. (but unknown to us).\n',r.ForValidation.UnknownBias);
disp('---------------------------------------')
end
%-------------------------------------------------------------
function B = EstimateBiases(r)
% gyroscopes data.
ggm = r.GyrosNoisyMesurements;  % Raw measurements from gyros.
ttg = r.timestamps;             % sampling times of gyros.

% sporadic attitude measurements.
AAm=r.MeasuredAttitudes;          % sporadic measurements of attitude (in degrees).

% AAm(1,:) :  roll measurements.
% AAm(2,:) :  pitch measurements.
% AAm(3,:) :  yaw measurements.
% (in Global CF)

tta=r.TimesMeasuredAttitudes;    %times at which the attitude was measured.
% AAm(:,k) was measured at time = tta(k).

iia=r.sampleIndexesA;             
% indexes of discrete times at which we got those sporadic attitude measurements. 
% tta = ttg(iia);  %<------- attitudes were measured at certain IMU sampling times.
% you decide if you want to use those indexes or not. 

A0 = r.InitialAttitude;          % initial attitude (at time 0).

PlotGyrosMeasurements(ggm,ttg);
PlotAttitudeSamples(tta,AAm,A0);

B=MyEstimateBiases(ggm,ttg,AAm,tta,iia,A0);   % your implementation.

end
%.................................................

% some plots to show the data.
function PlotGyrosMeasurements(ggm,ttg)
    figure(10) ; clf();
    subplot(211);
    plot(ttg,ggm(1,:),'r'); hold on; plot(ttg,ggm(2,:),'g'); plot(ttg,ggm(3,:),'b'); 
    title('measurements from gyros (noisy) (deg/second)');
    xlabel('time (seconds)');
    ylabel('measurements (deg/second)');
    legend({'gx','gy','gz'});
    pause(0.01);
end
%.................................................
function PlotAttitudeSamples(tta,aa,a0)

subplot(212);
plot(tta,aa(1,:),'*r'); hold on;plot(tta,aa(2,:),'*g');plot(tta,aa(3,:),'*b');
%plot(0,a0(1),'or'); hold on;plot(0,a0(2),'og');plot(0,a0(3),'ob');
xlabel('time (seconds)');
ylabel('(degrees)');
legend({'roll','pitch','yaw'});
title('few measured attitudes');
end
%.................................................

function B=MyEstimateBiases(ggm,ttg,AAm,tta,iia,A0);   % your implementation.
    AAm_rad = deg2rad(AAm);
    A0_rad = deg2rad(A0);
    
    % Convert gyro measurements to rad/s
    ggm_rad = deg2rad(ggm);
    
    % Set optimization options
    options = optimset('Display', 'iter', 'MaxIter', 100, 'TolFun', 1e-6);
    
    % Initial guess for biases (zero)
    initial_guess = zeros(3, 1);
    
    % Start timer
    tic;
    
    % Run optimization to find biases that minimize attitude error
    B = fminsearch(@(b) attitudeErrorCost(b, ggm_rad, ttg, AAm_rad, tta, A0_rad),initial_guess, options);
    
    % Calculate processing time
    processing_time = toc;
    
    % Convert biases back to deg/s for display
    B_deg = rad2deg(B);
    
    % Display results
    disp('---------------------------------------');
    fprintf('Estimated biases: [%.2f, %.2f, %.2f] deg/s\n', B_deg(1), B_deg(2), B_deg(3));
    fprintf('Processing time: %.2f seconds\n', processing_time);
    disp('---------------------------------------');pause(0.5);
  % to be implemented by you.
end
function cost = attitudeErrorCost(bias, gyro_meas, gyro_times, att_meas, att_times, att_init)
    % Remove bias from gyro measurements
    corrected_gyro = gyro_meas - bias;
    
    % Initialize attitude
    current_att = att_init;
    
    % Initialize cost
    cost = 0;
    
    % Index for next attitude measurement
    next_att_idx = 1;
    num_att_meas = size(att_meas, 2);
    
    % Integrate gyro measurements through time
    for i = 2:length(gyro_times)
        % Time step
        dt = gyro_times(i) - gyro_times(i-1);
        
        % Integrate one step
        current_att = integrateAttitude(corrected_gyro(:,i), dt, current_att);
        
        % Check if we've reached an attitude measurement time
        if next_att_idx <= num_att_meas && gyro_times(i) >= att_times(next_att_idx)
            
            % Calculate error between estimated and measured attitude
            error = current_att - att_meas(:,next_att_idx);
            
            % Add squared error to cost (weight all components equally)
            cost = cost + sum(error.^2);
            
            next_att_idx = next_att_idx + 1;
        end
    end
end

function new_att = integrateAttitude(gyro_rates, dt, current_att)
    % Simple attitude integration using Euler method
    % current_att = [roll; pitch; yaw] in radians
    
    roll = current_att(1);
    pitch = current_att(2);
    
    % Attitude dynamics equations
    roll_rate = gyro_rates(1) + sin(roll)*tan(pitch)*gyro_rates(2) + cos(roll)*tan(pitch)*gyro_rates(3);
    pitch_rate = cos(roll)*gyro_rates(2) - sin(roll)*gyro_rates(3);
    yaw_rate = sin(roll)/cos(pitch)*gyro_rates(2) + cos(roll)/cos(pitch)*gyro_rates(3);
    
    % Euler integration
    new_att = current_att + dt * [roll_rate; pitch_rate; yaw_rate];
end