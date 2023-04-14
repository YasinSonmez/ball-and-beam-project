classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html

        % Physical/simulation parameters
        g = 9.81;
        r_arm = 0.0254;
        L = 0.4255;
        K = 1.5;
        tau = 0.025;
        theta_saturation = 56 * pi / 180;

        % Observer parameters
        A = [0 1; 0 0];
        F = [0 1; 0 0];
        G = [0; 1];
        C = [1 0];
        H = [1 0];
        K1 = [2; 1];
        K2 = [2; 1];

        ddtheta_prev = 0;
        dtheta_prev = 0;
        % These variables may change every step
        t_prev = -1;
        u = 0;
        theta_d = 0;
        a_ball_ref_prev = 0;
        j_ball_ref_prev = 0;
        theta_ref_prev = 0;
        dtheta_ref_prev = 0;
        x_obs = [0 0 0 0]';
        V = 0;

        % To be initialized
        alpha
        beta

    end
    methods(Access = protected)
        function setupImpl(obj)
           obj.alpha = (5/7)*obj.g*obj.r_arm/obj.L;
           obj.beta = 5*obj.r_arm^2/7/obj.L^2;
        end

        function K = lqr_bb_model_2(obj, A,B,Q,R)
            R_inv = inv(R);
            G = B * R_inv * B';
        
            z_11 = A;
            z_12 = -G;
            z_21 = - Q;
            z_22 = -A';
        
            z = [z_11 z_12; z_21 z_22];
            w = z;
            w_prev = w;
            for i = 1:1000
                w = w - 0.5*(w-inv(w));
                %disp(norm(w-w_prev,"fro"));
            end 
        
           %n = size(A);
           n=4;
           
            w_11 = w(1:n, 1:n);
            w_12 = w(1:n, n+1:2*n);
            w_21 = w(n+1:2*n, 1:n);
            w_22 = w(n+1:2*n, n+1:2*n);
        
            M = [w_12;w_22 + eye(n)];
            N = [w_11 + eye(n); w_21];
        
            P = M\-N;
            K = R_inv*B'*P;
        end

        function V_servo = stepImpl(obj, t, p_ball, theta)
        % This is the main function called every iteration. You have to implement
        % the controller in this function, bu you are not allowed to
        % change the signature of this function. 
        % Input arguments:
        %   t: current time
        %   p_ball: position of the ball provided by the ball position sensor (m)
        %
        %   theta: servo motor angle provided by the encoder of the motor (rad)
        % Output:
        %   V_servo: voltage to the servo input.        
            %% Sample Controller: Simple Proportional Controller

            if t == obj.t_prev
                V_servo = obj.V;
                return
            end

            % Time step
            delta_t = t-obj.t_prev;

            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            
            % Change of variables
            p_ball_ref = p_ball_ref-obj.L/2;
            p_ball = p_ball-obj.L/2;
            
            % Extract x_2 and x_4 using observer
            w = obj.x_obs(1:2);
            z = obj.x_obs(3:4);
            phi = [0; obj.alpha*sin(theta)+obj.beta*p_ball*z(2)^2*cos(theta)^2];
            w(1) = w(1)-obj.L/2;
            w = w + (obj.A*w+phi+obj.K1*(p_ball-obj.C*w))*delta_t;
            w(1) = w(1)+obj.L/2;
            z = z + (obj.F*z+obj.G*obj.u+obj.K2*(theta-obj.H*z))*delta_t;
            obj.x_obs = [w; z];
            v_ball = w(2);
            dtheta = z(2);

            ddtheta = (dtheta - obj.dtheta_prev)/delta_t;
            
            theta_ref = obj.theta_ref_prev;
            dtheta_ref = obj.dtheta_ref_prev;
            theta_ref = 0;
%             theta_ref = theta_ref ...
%                 + sign(dtheta_ref)*sqrt((-a_ball_ref+obj.alpha*sin(theta_ref))/(obj.beta*(-p_ball_ref)*cos(theta_ref)^2))*delta_t;
            %theta_ref = min(theta_ref, theta_saturation);
            %theta_ref = max(theta_ref, -theta_saturation);
%             dtheta_ref = sign(dtheta_ref)*sqrt((-a_ball_ref+obj.alpha*sin(theta_ref))/(obj.beta*(-p_ball_ref)*cos(theta_ref)^2));
            dtheta_ref = 0;
            u_ref = (obj.tau*(dtheta_ref-obj.dtheta_ref_prev)/delta_t + dtheta_ref)/obj.K;
            
            x_ref = [p_ball_ref; v_ball_ref; theta_ref; dtheta_ref];
            x_states = [p_ball; v_ball; theta; dtheta];
            
            A_lqr = [0 1 0 0; ...
                obj.beta*dtheta_ref^2*cos(theta_ref)^2 0 ...
                obj.alpha*cos(theta_ref)+2*obj.beta*p_ball_ref*dtheta_ref^2*cos(theta_ref)*sin(theta_ref) 2*obj.beta*p_ball_ref*dtheta_ref^2*cos(theta_ref)^2; ...
                0 0 0 1; ...
                0 0 0 -1/obj.tau];
            B_lqr = [0; 0; 0; obj.K/obj.tau];

            if p_ball < -0.30
                Q = [18000 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
                R = 1;
            elseif p_ball > -0.11
                Q = [18000 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
                R = 1;
            else
                Q = [18000 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
                R = 5;
           end

            if a_ball_ref == 0 %Q and R cost for square wave
                Q = [1800 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
                R = 5; 
            end 
            
            K_lqr = lqr_bb_model_2(obj, A_lqr, B_lqr, Q, R);
            
            
            if theta > obj.theta_saturation
                V_servo = ddtheta*obj.tau/obj.K;
            elseif theta < -obj.theta_saturation
                V_servo = ddtheta*obj.tau/obj.K;
            else 
                V_servo = -K_lqr*(x_states - x_ref) + u_ref;
            end

            obj.u = -dtheta/obj.tau+V_servo*obj.K/obj.tau;

            % Update variables
            obj.dtheta_prev = dtheta;
            obj.ddtheta_prev = ddtheta;
            obj.t_prev = t;
            obj.theta_ref_prev = theta_ref;
            obj.dtheta_ref_prev = dtheta_ref;
            obj.V = V_servo;
        end
    end
    
    methods(Access = public)
        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d, x_obs] = stepController(obj, t, p_ball, theta)
            if t == 0
                setupImpl(obj)
            end
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
            x_obs = obj.x_obs;
        end
    end
    
end
