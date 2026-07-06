function [outputArg1,outputArg2] = convert_to_unstructured(varargin,nx,ny,nz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nx = Nx;
ny = Ny;
nz = Nz;
n = nx*ny*nz;

indX = [1:nx]';
indY = [1:ny]';
indZ = [1:nz]';


indXY = [repelem(indX, size(indY,1),1), repmat(indY, size(indX,1),1)];
indXYZ = [repmat(indXY, size(indZ,1),1), repelem(indZ, size(indXY,1),1)];
linear_ind = sub2ind([nx, ny, nz], indXYZ(:,1),indXYZ(:,2),indXYZ(:,3));
C = reshape(linear_ind,nx,ny,nz);
back = reshape(C,[],1,1);
Fx = [reshape(C(1:nx-1,:,:),[],1,1) reshape(C(2:nx,:,:),[],1,1)];
Fy = [reshape(C(:,1:ny-1,:),[],1,1) reshape(C(:,2:ny,:),[],1,1)];
Fz = [reshape(C(:,:,1:nz-1),[],1,1) reshape(C(:,:,2:nz),[],1,1)];

%% preparing face ind

    % Ugx = zeros(Nx+1, Ny, Nz);
    %                 Ugy = zeros(Nx, Ny+1, Nz);
    %                 Ugz = zeros(Nx, Ny, Nz+1);


fIndX = zeros(nx+1, ny, nz);
fIndX(2:nx,:) = reshape([1:size(Fx,1)]',nx-1, 1, []);
fIndX_flat = reshape(fIndX,[],1,1);

F = [Fx;Fy;Fz];


                 % U.x(2:Nx,:,:) = Grid.Tx(2:Nx,:,:) .* ( (P(1:Nx-1,:,:)-P(2:Nx,:,:)) - Ugx(2:Nx,:,:) );
                 %    U.y(:,2:Ny,:) = Grid.Ty(:,2:Ny,:) .* ( (P(:,1:Ny-1,:)-P(:,2:Ny,:)) - Ugy(:,2:Ny,:) );
                 %    U.z(:,:,2:Nz) = Grid.Tz(:,:,2:Nz) .* ( (P(:,:,1:Nz-1)-P(:,:,2:Nz)) - Ugz(:,:,2:Nz) );


if size(varargin,1) == 1
    if isstruct(varargin)
        arg = varargin;
    else
        arg = varargin{1};
    end
        xx = arg.x;
        yy = arg.y;
        zz = arg.z;
        ff = permute(zz,[1 3 2]);        
        ff = ff(:,2:end);
        shape_shift(zz, nx,ny,nz)
else
        arg = varargin;
        shape_shift(arg, nx,ny,nz)
end


end
function output = shape_shift(arg, nx,ny,nz)
    num_col = size(arg,2);
    n = nz*ny*nx;
    output = zeros(n,num_col);
    for index = 1:n
        output(:,index) = reshape(arg(:,index),[],1,1)
    end
        

end