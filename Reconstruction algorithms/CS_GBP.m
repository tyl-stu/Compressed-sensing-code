function hat_x=CS_GBP(y,T_Mat,m)
% y=T_Mat*x, T_Mat is n-by-m
% y - measurements
% T_Mat - combination of random matrix and sparse representation basis,A
% m - size of the original signal
% the sparsity is length(y)/4

cnt=length(y); 

x=y';                         % x is row vector
D=[T_Mat,-T_Mat]';            % each row of D is a atom


OPT_VERBOSE_ON = 0; % commentary level
epsilon = 1.0e-6; % error tolerance
[n, d] = size(D); % problem size
pinv_epsilon = 1.0e-12; %% epsilon for PINV_ITER
exit_flag = 0; % default exit_flag

% Compute the initial atom
t = 1; % iteration index
[alpha_0, k] = max(D*x');
i_list = [k];
alpha_list = [alpha_0];
n_atoms = length(i_list); % number of atoms in the current representation
normal = x / sqrt(sum(x.^2)); % the normal to the hyperplane
xApprox = alpha_0*D(k,:); % the current approximation to x
d_H = dot(D(k,:), normal); % the distance of the hyperplane to the origin
xApprox_H = (d_H / dot(xApprox, normal))*xApprox; % xApprox projected onto the hyperplane
err_vec = x - xApprox; %
v = err_vec - dot(err_vec, normal)*normal;
v = v / sqrt(sum(v.^2));
%
Psi = [D(k,:)];
Psi_biortho = Psi;

% Iterate
for times=1:cnt/4
    % Exit when error is acceptable
    error = sqrt(sum(err_vec.^2));
    if OPT_VERBOSE_ON,
        fprintf('t = %i; i_t = %i; nAtoms = %i; error = %f\n', t, k, n_atoms, error);
    end

    if error < epsilon,
        fprintf('... GBP complete.\n');
        break;
    elseif isnan(error),
        fprintf('GBP -- Error: Exiting.\n');
        exit_flag = -2;
    end

    % Exit when number of atoms equals dimension
    if n_atoms==d,
        fprintf('Basis found.\n');
        break;
    end
    % Otherwise, next iteration
    t = t + 1;

    % Project atoms onto n-v plane
    indices = setdiff(1:n, i_list);  %delete the elements of i_list from 1:n
    D_n0 = D(indices,:) * normal';  % select the next atom from the residuals
    D_v0 = D(indices,:) * v';
    D_n = D_n0 - repmat(dot(xApprox_H, normal), length(indices), 1);
    D_v = D_v0 - repmat(dot(xApprox_H, v), length(indices), 1);
    % Compute angle
    D_theta = atan2(-D_n, D_v);
    D_theta = D_theta + (D_theta < 0)*2*pi; % re-center for min()

    % Order by angle and select the next intersected atom
    [junk, i_Dsorted] = sort(D_theta);
    k = indices(i_Dsorted(1));
    psi_k = D(k,:);

    % Compute auxiliary variables and coefficients
    % Compute orthogonal component of psi_k
    psi_k_ortho = psi_k;
    for i = 1:n_atoms,
        beta(i) = dot(psi_k, Psi_biortho(i,:));
        psi_k_ortho = psi_k_ortho - beta(i)*Psi(i,:);
    end
    % Compute new biorthogonal vector and adjust existing biorthogonal vectors
    psi_k_biortho = psi_k_ortho / sum(psi_k_ortho.^2);
    for i = 1:n_atoms,
        Psi_biortho(i,:) = Psi_biortho(i,:) - beta(i)*psi_k_biortho;
    end
    % Compute new coefficient and adjust existing coefficients
    alpha_k = dot(x, psi_k_biortho);
    for i = 1:n_atoms,
        alpha_list(i) = alpha_list(i) - beta(i)*alpha_k;
    end

    % Check if the new representation is convex (with respect to the new coefficient)
    % If not, there is a numerical error.
    if alpha_k < 0,
        fprintf('GBP -- Error: Hyperplane rotation failed. Exiting. \n');
        exit_flag = -1;
        return;
    end

    % Update lists
    Psi = [Psi; psi_k];
   
    Psi_biortho = [Psi_biortho; psi_k_biortho];
    alpha_list = [alpha_list, alpha_k];
    i_list = [i_list, k];
    n_atoms = n_atoms + 1;

    %Test biorthogonality
    AAplus = Psi*Psi_biortho';
    Ierr = max(max(abs(AAplus - eye(size(AAplus)))));
    if Ierr > epsilon,
        % run pinv_iter with epsilon=epsilon and c=1
        X = Psi_biortho';
        while Ierr > epsilon,
            X = X * (2*eye(size(AAplus)) - AAplus);
            AAplus = Psi*X;
            Ierr = max(max(abs(AAplus - eye(size(AAplus)))));
            %fprintf('... pinv: I error: %f\n', Ierr);
        end
        % pinv_iter complete
        Psi_biortho = X';
        alpha_list = (Psi_biortho*x')';
    end

    % Test if conv(Psi) contains extraneous atoms,
    % otherwise, continue wrapping
    while sum(alpha_list<=0)~=0,
        if OPT_VERBOSE_ON,
            fprintf('... Deleting atom. \n');
        end
        % Find an extraneous atom
        [alpha_j, j] = min(alpha_list);
        % Delete it
        psi_j_biortho = Psi_biortho(j,:);
        gamma = Psi_biortho*psi_j_biortho' / sum(psi_j_biortho.^2);
        Psi = Psi([1:j-1,j+1:n_atoms],:);
        
        Psi_biortho = Psi_biortho - gamma*psi_j_biortho;
        Psi_biortho = Psi_biortho([1:j-1,j+1:n_atoms],:);

        % Test biorthogonality
        AAplus = Psi*Psi_biortho';
        Ierr = max(max(abs(AAplus - eye(size(AAplus)))));
        if Ierr > epsilon,
            % run pinv_iter with epsilon=epsilon and c=1
            X = Psi_biortho';
            while Ierr > epsilon,
                X = X * (2*eye(size(AAplus)) - AAplus);
                AAplus = Psi*X;
                Ierr = max(max(abs(AAplus - eye(size(AAplus)))));
                %fprintf('... pinv: I error: %f\n', Ierr);
            end
            % pinv_iter complete
            Psi_biortho = X';
            alpha_list = (Psi_biortho*x')';
        else
            alpha_list = alpha_list - alpha_j*gamma';
            alpha_list = alpha_list([1:j-1,j+1:n_atoms]);
        end

        % Update the representation
        i_list = i_list([1:j-1,j+1:n_atoms]);
        n_atoms = length(i_list);
    end

    % Compute the new approximation
    xApprox = alpha_list * Psi;
    % Next hyperplane
    normal = (-D_n(i_Dsorted(1)))*v + (D_v(i_Dsorted(1)))*normal;
    normal = normal / sqrt(sum(normal.^2));
    d_H = dot(D(i_list(1),:), normal);
    xApprox_H = (d_H / dot(xApprox, normal))*xApprox;
    err_vec = x - xApprox;
    v = err_vec - dot(err_vec, normal)*normal; 
    v = v / sqrt(sum(v.^2));
end

hat_x=zeros(1,size(D,1));
hat_x(i_list)=alpha_list;         % recovered sparse representation
end