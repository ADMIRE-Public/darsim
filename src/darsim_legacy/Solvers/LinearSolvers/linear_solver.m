% Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Last modified: 21 May 2025 by Mathias Erdtmann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver < handle
    properties
        Name
        Tol
        Maxit
        Iter
        CondNumber
        %Preconditioner
    end
    methods
        function obj = linear_solver(name, tol, maxit)
            obj.Name = name;
            obj.Tol = tol;
            obj.Maxit = maxit;
        end
        function SetUp(obj, ProductionSystem, DiscretizationModel, Residual)
            %obj.Preconditioner.Setup(ProductionSystem, DiscretizationModel);
        end
        function x = Solve(obj, A, rhs)
            %obj.CondNumber = condest(A);
            obj.CondNumber = 1;
            start = tic;       
            switch (obj.Name)
                case('gmres')
                    % Set-up ilu preconditioner
                    setup.type = 'nofill';
                    setup.milu = 'off';
                    setup.droptol = 0.1;
                    Diagonal = diag(A);
                    i=find(Diagonal==0);
                    if ~isempty(i)
                        disp(i);
                    end
                    [L, U] = ilu(A, setup);
                    [x, flag, relres, obj.Iter] = gmres(A, rhs, [], obj.Tol, obj.Maxit, L, U);
                case('bicg')
                    % Set-up ilu preconditioner
                    setup.type = 'nofill';
                    setup.milu = 'off';
                    setup.droptol = 0.1;
                    Diagonal = diag(A);
                    i=find(Diagonal==0);
                    if ~isempty(i)
                        disp(i);
                    end
                    [L, U] = ilu(A, setup);
                    [x, flag, relres, obj.Iter] = bicg(A, rhs, obj.Tol, obj.Maxit, L, U);
                case('direct')
                    flag = 0;
                    if isempty(gcp('nocreate')) % find if the multicore pool was initiated
                        x = A\rhs;
                    % Parallel linear solver
                    % Written by Mathias Erdtmann
                    else % Parallel pool is active with >1 workers.
                        spmd
                            localRows = codistributed(A); % Distribute rows across workers
                            localB = codistributed(rhs);    % Distribute identity matrix
		    	            localInverse = localRows \ localB; % Solve for local part of the inverse
                        end
		                x = gather(localInverse);
                    end
                case('symmlq')
                    [x, flag, relres, obj.Iter] = symmlq(A'*A,A'*rhs,obj.Tol,obj.Maxit);
                otherwise
                    error('unsupported linear solver type');
            end
            % disp(['LS time:', num2str(toc(start))]);
            if flag ~= 0
                disp(['WARNING: Linear solver did not converge. The residual norm is ', num2str(relres)]);
            end
        end
    end
end