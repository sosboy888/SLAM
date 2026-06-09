classdef PlatformPredictionEdge < g2o.core.BaseBinaryEdge
    % PlatformPredictionEdge summary of PlatformPredictionEdge
    %
    % This class stores the factor representing the process model which
    % transforms the state from timestep k to k+1
    %
    % The process model is as follows.
    %
    % Define the rotation vector
    %
    %   M = dT * [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0;0 0 1];
    %
    % The new state is predicted from 
    %
    %   x_(k+1) = x_(k) + M * [vx;vy;theta]
    %
    % Note in this case the measurement is actually the mean of the process
    % noise. It has a value of 0. The error vector is given by
    %
    % e(x,z) = inv(M) * (x_(k+1) - x_(k))
    %
    % Note this requires estimates from two vertices - x_(k) and x_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k) and slot 2 contains x_(k+1).
    
    properties(Access = protected)
        % The length of the time step
        dT;
    end
    
    methods(Access = public)
        function obj = PlatformPredictionEdge(dT)
            % PlatformPredictionEdge for PlatformPredictionEdge
            %
            % Syntax:
            %   obj = PlatformPredictionEdge(dT);
            %
            % Description:
            %   Creates an instance of the PlatformPredictionEdge object.
            %   This predicts the state from one timestep to the next. The
            %   length of the prediction interval is dT.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a PlatformPredictionEdge

            assert(dT >= 0);
            obj = obj@g2o.core.BaseBinaryEdge(3);            
            obj.dT = dT;
        end
       
        function initialEstimate(obj)
            % INITIALESTIMATE Compute the initial estimate of a platform.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial estimate of the platform x_(k+1) given
            %   an estimate of the platform at time x_(k) and the control
            %   input u_(k+1)

            xk = obj.edgeVertices{1}.x;
            theta = xk(3);
        
            dT = obj.dT;
            if dT < 1e-12
                % No time has passed, best we can do is copy forward
                obj.edgeVertices{2}.x = xk;
                return;
            end
        
            % Control input/odometry is the measurement z
            u = obj.z; %3x1
        
            M = dT * [ cos(theta) -sin(theta) 0;
                       sin(theta)  cos(theta) 0;
                       0           0          1];
        
            obj.edgeVertices{2}.x = xk + M * u;
        end
        
        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the measurement and the parameter state in the
            %   vertex. Note the error enters in a nonlinear manner, so the
            %   equation has to be rearranged to make the error the subject
            %   of the formulat
                       
            xk   = obj.edgeVertices{1}.x;
            xkp1 = obj.edgeVertices{2}.x;
        
            theta = xk(3);
            dT = obj.dT;
            if dT < 1e-12
                obj.errorZ = zeros(3,1);
                return;
            end
        
            u = obj.z; %3x1
        
            % Use the closed form inverse
            c = cos(theta); s = sin(theta);
            Minv = (1/dT) * [ c  s  0;
                             -s  c  0;
                              0  0  1];
        
            obj.errorZ = Minv * (xkp1 - xk) - u;
        end
        
        % Compute the Jacobians
        function linearizeOplus(obj)
            % LINEARIZEOPLUS Compute the Jacobians for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the Jacobians for the edge. Since we have two
            %   vertices which contribute to the edge, the Jacobians with
            %   respect to both of them must be computed.
            %

            xk   = obj.edgeVertices{1}.x;
            xkp1 = obj.edgeVertices{2}.x;
        
            theta = xk(3);
            dT = obj.dT;
            if dT < 1e-12
                obj.J{1} = zeros(3,3);
                obj.J{2} = zeros(3,3);
                return;
            end
        
            c = cos(theta); s = sin(theta);
        
            Minv = (1/dT) * [ c  s  0;
                             -s  c  0;
                              0  0  1];
        
            % Base parts(from Minv*(xkp1-xk))
            J1 = -Minv;
            J2 =  Minv;
        
            % Minv depends on theta in xk
            dx = xkp1(1) - xk(1);
            dy = xkp1(2) - xk(2);
        
            dMinv_dtheta = (1/dT) * [ -s  c  0;
                                     -c -s  0;
                                      0  0  0];
        
            dtheta_term = dMinv_dtheta * [dx; dy; (xkp1(3)-xk(3))];
        
            % Only the 3rd column(theta) gets this extra term
            J1(:,3) = J1(:,3) + dtheta_term;
        
            obj.J{1} = J1;
            obj.J{2} = J2;
        end
    end    
end