function [Xi] = SS_Soln_Aug(obj,~,As,Bs,ts,u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



if nargin==1 
    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts;
    u = obj.u;
    keep_SS = false;
elseif nargin == 2
    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts;
    u = obj.u;
elseif nargin ~= 6
    fprintf('There were not enough or too many inputs provided.')
end

ti= ts;

depends = obj.Converter.Topology.Parser.dependsAB;


    ns = size(As,2);
    pi = size(Bs,2);
    nsub = length(ti);
%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This sets up dummy dependent states that will be added back in
    % and be incorrect but they will be fixed by CorrectXs at the end
    % of this function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Is = [eye(size(As,1)-size(depends,1)) zeros(size(As,1)-size(depends,1),size(depends,1)) ; eye(size(depends)) zeros(size(depends,1))  ];



  %% Set up Augmented State Space
    Phi = eye(ns+1); %State Transition Matrix; multiplied out in next loop
    for i = 1:nsub
        Aas(:,:,i) = [As(:,:,i), Bs(:,:,i)*u; 
                        zeros(1,ns+1)];
        Phi = expm(Aas(:,:,i)*ti(i))*Phi;
    end
    
    %% Get rid of Dependent States
    % (removing dependnet states improves numerical accuracy)
    depStates = find(sum(abs(eye(ns)-Is(:,:,1)),2));
    nd = length(depStates);
    Phi(depStates, :, :) = [];
    Phi(:,depStates,:) = [];
    
    %% Solve Augmented State Space without dependent states
    x0b = null(eye(ns-nd+1) - Phi);
    x0b  = x0b(:,end);
    X0 = x0b/x0b(end);
    
    %% Put dependent states back in
    % Putting them back in makes indexing to find states of interest easier, later
    for i = 1:length(depStates)
       X0 = [X0(1:depStates(i)-1); 
             0; 
             X0(depStates(i):end)];
    end
    X0 = [Is(:,:,1)*X0(1:end-1); X0(end)];
    

    %% Integrate all states and outputs during each subinterval
    Xi = zeros(2*ns+1,nsub+1);
    Xi(:,1) = [X0; zeros(ns,1)];
    for i = 1:nsub
        % Triple-augmented matrix; in addition to the eye(ns) term that
        % integrates the states, this adds in the appropriate C and D
        % matrices such that the last mo entries in the new state vector
        % are the integrals of each of the outputs
        Aais(:,:,i) = [Aas(:,:,i), zeros(ns+1, ns);
                        eye(ns), zeros(ns, ns+1)];
        Xi(:,i+1) = expm(Aais(:,:,i)*ti(i))*Xi(:,i);
    end
    
    Xi(size(As,1)+1:end,:) = [];
    obj.Converter.Topology.Parser.StateVarIndex();
    [Xi] = obj.CorrectXs(0,Xi);
    
    
    
    %% Check that solution is valid
    assert(norm(Xi(1:ns-nd,1)-Xi(1:ns-nd,end)) < 1e-6)
 

    
    

end