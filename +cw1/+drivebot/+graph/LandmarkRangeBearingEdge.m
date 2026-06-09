classdef LandmarkRangeBearingEdge < g2o.core.BaseBinaryEdge
    % LandmarkRangeBearingEdge summary of LandmarkRangeBearingEdge
    %
    % This class stores an edge which represents the factor for observing
    % the range and bearing of a landmark from the vehicle. Note that the
    % sensor is fixed to the platform.
    %
    % The measurement model is
    %
    %    z_(k+1)=h[x_(k+1)]+w_(k+1)
    %
    % The measurements are r_(k+1) and beta_(k+1) and are given as follows.
    % The sensor is at (lx, ly).
    %
    %    dx = lx - x_(k+1); dy = ly - y_(k+1)
    %
    %    r(k+1) = sqrt(dx^2+dy^2)
    %    beta(k+1) = atan2(dy, dx) - theta_(k+1)
    %
    % The error term
    %    e(x,z) = z(k+1) - h[x(k+1)]
    %
    % However, remember that angle wrapping is required, so you will need
    % to handle this appropriately in compute error.
    %
    % Note this requires estimates from two vertices - x_(k+1) and l_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k+1) and slot 2 contains l_(k+1).
    
    methods(Access = public)
    
        function obj = LandmarkRangeBearingEdge()
            % LandmarkRangeBearingEdge for LandmarkRangeBearingEdge
            %
            % Syntax:
            %   obj = LandmarkRangeBearingEdge();
            %
            % Description:
            %   Creates an instance of the LandmarkRangeBearingEdge object.
            %   Note we feed in to the constructor the landmark position.
            %   This is to show there is another way to implement this
            %   functionality from the range bearing edge from activity 3.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a LandmarkRangeBearingEdge.

            obj = obj@g2o.core.BaseBinaryEdge(2);
        end
        
        function initialEstimate(obj)
            % INITIALESTIMATE Compute the initial estimate of the landmark.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial estimate of the landmark given the
            %   platform pose and observation.

            pose = obj.edgeVertices{1}.x;
            landmark = obj.edgeVertices{2}.x;
            
            % Only initialise if landmark is not set
            if any(isnan(landmark))
            
                r = obj.z(1);
                beta = obj.z(2);
            
                x = pose(1);
                y = pose(2);
                psi = pose(3);
            
                alpha = psi + beta;
            
                landmark = [
                    x + r*cos(alpha);
                    y + r*sin(alpha)];
            
                obj.edgeVertices{2}.x = landmark;
            
            end
        end
        
        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the predicted and actual range-bearing measurement.

            pose = obj.edgeVertices{1}.x;
            landmark = obj.edgeVertices{2}.x;
            
            x = pose(1);
            y = pose(2);
            psi = pose(3);
            
            xi = landmark(1);
            yi = landmark(2);
            
            dx = xi - x;
            dy = yi - y;
            
            range_pred = sqrt(dx^2 + dy^2);
            
            bearing_pred = atan2(dy, dx) - psi;
            
            error_range = obj.z(1) - range_pred;
            
            error_bearing = obj.z(2) - bearing_pred;
            
            % wrap angle
            error_bearing = mod(error_bearing + pi, 2*pi) - pi;
            
            obj.errorZ = [
                error_range;
                error_bearing];
        end
        
        function linearizeOplus(obj)
            % linearizeOplus Compute the Jacobian of the error in the edge.
            %
            % Syntax:
            %   obj.linearizeOplus();
            %
            % Description:
            %   Compute the Jacobian of the error function with respect to
            %   the vertex.
            %

            pose = obj.edgeVertices{1}.x;
            landmark = obj.edgeVertices{2}.x;
            
            x = pose(1);
            y = pose(2);
            
            xi = landmark(1);
            yi = landmark(2);
            
            dx = xi - x;
            dy = yi - y;
            
            q = dx^2 + dy^2;
            
            range = sqrt(q);
            
            % Jacobian wrt pose
            
            J_pose = [
            
             dx/range    dy/range     0;
            
             -dy/q       dx/q         1
            
            ];
            
            % Jacobian wrt landmark
            
            J_landmark = [
            
             -dx/range   -dy/range;
            
              dy/q       -dx/q
            
            ];
            
            obj.J{1} = J_pose;
            obj.J{2} = J_landmark;
        end        
    end
end
