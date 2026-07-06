function [intersectCoord, areCoplanar, indNeighbors] = Obtain_Tetrahedron_Polygon_Intersection(Plane_fracCell, TetrahedralGridData, ...
                                                                                               CellInnerFaces, CellBoundFaces, fintpoint, ...
                                                                                               fboupoint, monjuscell, Im, Epsilon)
%
% Os planos estão definidos como: -----------------------------------------
% pi_1: N1(:,1)*x + N1(:,2)*y + N1(:,3)*z = ( N1(:,1)*J(:,1) + N1(:,2)*J(:,2) + N1(:,3)*J(:,3) ) 
% pi_2: N2(:,1)*x + N2(:,2)*y + N2(:,3)*z = ( N2(:,1)*J(:,1) + N2(:,2)*J(:,2) + N2(:,3)*J(:,3) )
% pi_3: N3(:,1)*x + N3(:,2)*y + N3(:,3)*z = ( N3(:,1)*J(:,1) + N3(:,2)*J(:,2) + N3(:,3)*J(:,3) )
% pi_4: N4(:,1)*x + N4(:,2)*y + N4(:,3)*z = ( N4(:,1)*I(:,1) + N4(:,2)*I(:,2) + N4(:,3)*I(:,3) ) 
% pi_f: Nf(:,1)*x + Nf(:,2)*y + Nf(:,3)*z = ( Nf(:,1)*Jf(:,1) + Nf(:,2)*Jf(:,2) + Nf(:,3)*Jf(:,3) )
%--------------------------------------------------------------------------
E = -log10(Epsilon);

TetraCoord = TetrahedralGridData.Vertex.Coordinates( TetrahedralGridData.Cell.Vertices(Im,:) , : );
FractCoord = [ Plane_fracCell.PointA; Plane_fracCell.PointB; Plane_fracCell.PointC; Plane_fracCell.PointD ];

if Im==1, ini = 1; else ini = fintpoint(Im-1) + 1; end
TetraIFace = CellInnerFaces( ini : fintpoint(Im), 2 );
% TetraIFace = CellInnerFaces(CellInnerFaces(:,1)==Im,2);

if isempty(TetraIFace)==0
    TetraIFaceVert = TetrahedralGridData.Face.Inner.Vertices(TetraIFace,:);
    TetraIFaceNorm = TetrahedralGridData.Face.Inner.Normals(TetraIFace,:);
end

if Im==1, ini = 1; else ini = fboupoint(Im-1) + 1; end
TetraBFace = CellBoundFaces( ini : fboupoint(Im), 2 );
% TetraBFace = CellBoundFaces(CellBoundFaces(:,1)==Im,2);
if isempty(TetraBFace)==0
    TetraBFaceVert = TetrahedralGridData.Face.Bound.Vertices(TetraBFace,:);
    TetraBFaceNorm = TetrahedralGridData.Face.Bound.Normals(TetraBFace,:);
else
    TetraBFaceVert = [];
    TetraBFaceNorm = [];
end

TetraFaceVert = [ TetraIFaceVert; TetraBFaceVert ];
TetraFaceNorm = [ TetraIFaceNorm; TetraBFaceNorm ];

TetFacVertCoord1 = TetrahedralGridData.Vertex.Coordinates(TetraFaceVert(:,1),:);

FractNorm = cross(FractCoord(2,:) - FractCoord(1,:), FractCoord(3,:) - FractCoord(1,:)); FractNorm = FractNorm/norm(FractNorm);
FractNormMat = [ FractNorm; FractNorm; FractNorm; FractNorm ];

Frac2TetVecs = [ TetraCoord(1,:) - FractCoord(1,:); TetraCoord(2,:) - FractCoord(1,:); ...
       TetraCoord(3,:) - FractCoord(1,:); TetraCoord(4,:) - FractCoord(1,:) ];
   
veri1 = abs(sum(sign(round(sum(Frac2TetVecs.*FractNormMat,2),E))))<4;

if abs(sum(sign(round(sum(Frac2TetVecs.*FractNormMat,2),E))))==1
    
    areCoplanar = 1;
    
    vertface = TetrahedralGridData.Cell.Vertices(Im,find(sign(round(sum(Frac2TetVecs.*FractNormMat,2),E))==0));
    
    intersectCoord = TetrahedralGridData.Vertex.Coordinates( vertface, : );
    
    intersectCoord = intersectionpolys( intersectCoord, FractCoord );
    
elseif veri1==0
    
    intersectCoord = []; areCoplanar = 0;
    
else
    
    plan = [1 2; 1 3; 2 3; 2 4; 1 4; 4 3]; t = 1;

    for i=1:6

        M = [ TetraFaceNorm(plan(i,1),:); TetraFaceNorm(plan(i,2),:); FractNorm ];

        if ~((M(1,1)==M(2,1) && M(1,1)==M(3,1))||(M(1,2)==M(2,2) && M(1,2)==M(3,2))||(M(1,3)==M(2,3) && M(1,3)==M(3,3))) && rank(M,1e-8)==3

            B = [ TetraFaceNorm(plan(i,1),:)*TetFacVertCoord1(plan(i,1),:)'; TetraFaceNorm(plan(i,2),:)*TetFacVertCoord1(plan(i,2),:)'; FractNorm*FractCoord(1,:)'];

            intersectCoord(t,:) = (M\B)'; t = t + 1;

        end

    end
    
    intersectCoord = vertinside( intersectCoord, TetraFaceNorm, TetraFaceVert, TetrahedralGridData.Vertex.Coordinates, E );
    
    intersectCoord = intersectionpolys( intersectCoord, FractCoord );
    
    areCoplanar = 0;

end

todelem = monjuscell(TetraIFace,:);
outelpos = monjuscell(TetraIFace,:)~=Im;
indNeighbors = todelem(outelpos);

intersectCoord = unique(intersectCoord,'rows');

% if isempty(intersectCoord)==0
%     
%     plotintersec( TetraFaceVert, TetrahedralGridData.Vertex.Coordinates, FractCoord, intersectCoord );
% 
% end

end

function intersectCoord = vertinside( intersectCoord, TetraFaceNorm, TetraFaceVert, Coord, E )
%

for i=1:size(TetraFaceVert,1)
    
    MeiosdasFaces(i,:) = (1/3)*sum( Coord( TetraFaceVert(i,:), : ), 1 );
    
end

MeiodoTetra = 0.25*sum( MeiosdasFaces, 1 );

Meio2Faces = MeiosdasFaces - repmat( MeiodoTetra, [ size(MeiosdasFaces,1) 1 ] );

verisentnorm = sign(sum(Meio2Faces.*TetraFaceNorm,2));

TetraFaceNorm = TetraFaceNorm.*[ verisentnorm verisentnorm verisentnorm ];

for i=1:size( intersectCoord,1 )
    
    Meio2Int = repmat( intersectCoord(i,:), [ size(MeiosdasFaces,1) 1 ] ) - MeiosdasFaces;
    
    veriinside(i) = ~any(sign(round(sum(Meio2Int.*TetraFaceNorm,2),E))==1);
    
end

intersectCoord = intersectCoord(veriinside==1,:);

end

function intersectCoord = intersectionpolys( intersectCoord, FractCoord )
%

O = FractCoord(1,:); i = [ 1 0 0 ]; j = [ 0 1 0 ]; k = [ 0 0 1 ];

vecint = intersectCoord - O; vecfract = FractCoord - O;

lin = find( abs(vecfract(:,3))<1e-10 & ( abs(vecfract(:,1))>1e-10 | abs(vecfract(:,2))>1e-10 ) );

vecdir = vecfract(lin,:); vecdir = vecdir/norm(vecdir);

theta = round(atan(vecdir(2)/vecdir(1)),10);

R = [ cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1 ];

vecfractrot = vecfract*R; vecintrot = vecint*R;

nvfrot = cross( vecfractrot(2,:), vecfractrot(3,:) );

alfa = round(atan(nvfrot(3)/nvfrot(2)),10);

R = [ 1 0 0; 0 cos(alfa) -sin(alfa); 0 sin(alfa) cos(alfa) ];

vecfractrot = vecfractrot*R; vecintrot = vecintrot*R; vecintrot = reordint( vecintrot );
%--------------------------------------------------------------------------
vecintPlanCoord = [ ];
%--------------------------------------------------------------------------
xq = vecintrot(:,1); yq = vecintrot(:,3); xv = vecfractrot(:,1); yv = vecfractrot(:,3); 

[ inint, ~ ] = inpolygon( xq, yq, xv, yv );

if ~all(inint) && any(inint)
    
    xq = [ xq; xq(1) ]; yq = [ yq; yq(1) ];
    
    xv = [ xv; xv(1) ]; yv = [ yv; yv(1) ];
    
    [ xi, yi ] = polyxpoly( xq, yq, xv, yv );
    
    pontosdentro = [ xi zeros(size(xi)) yi ];
    
    vecintPlanCoord = [ vecintPlanCoord; pontosdentro ];
    
end

vecintPlanCoord = [ vecintPlanCoord; vecintrot(inint,:) ];
%--------------------------------------------------------------------------
xv = vecintrot(:,1); yv = vecintrot(:,3); xq = vecfractrot(:,1); yq = vecfractrot(:,3); 

[ infrac, ~ ] = inpolygon( xq, yq, xv, yv );

if ~all(infrac) && any(infrac)
    
    xq = [ xq; xq(1) ]; yq = [ yq; yq(1) ];
    
    xv = [ xv; xv(1) ]; yv = [ yv; yv(1) ];
    
    [ xi, yi ] = polyxpoly( xq, yq, xv, yv );
    
    pontosdentro = [ xi zeros(size(xi)) yi ];
    
    vecintPlanCoord = [ vecintPlanCoord; pontosdentro ];
    
end

vecintPlanCoord = [ vecintPlanCoord; vecfractrot(infrac,:) ];
%--------------------------------------------------------------------------
vecintPlanCoord = unique( vecintPlanCoord, 'rows' );

vecintPlanCoord = reordint( vecintPlanCoord );
%--------------------------------------------------------------------------

R = [ 1 0 0; 0 cos(-alfa) -sin(-alfa); 0 sin(-alfa) cos(-alfa) ];

vecintCoord = vecintPlanCoord*R;

R = [ cos(-theta) -sin(-theta) 0; sin(-theta) cos(-theta) 0; 0 0 1 ];

intersectCoord = vecintCoord*R + O;
%--------------------------------------------------------------------------

end

function plotintersec( TetraFaceVert, Coord, FractCoord, intersectCoord )
%

TetraFaceVertRep = [ TetraFaceVert TetraFaceVert(:,1) ];

for i=1:size(TetraFaceVert,1)
   
    TetraCoordRep = Coord(TetraFaceVertRep(i,:),:);
    
    plot3(TetraCoordRep(:,1),TetraCoordRep(:,2),TetraCoordRep(:,3),'-k')
    
    hold on
    
end

FractCoordRep = [ FractCoord; FractCoord(1,:) ];

plot3(FractCoordRep(:,1),FractCoordRep(:,2),FractCoordRep(:,3),'-r')

intersectCoord = [ intersectCoord; intersectCoord(1,:) ];

plot3(intersectCoord(:,1),intersectCoord(:,2),intersectCoord(:,3),'-b')

hold off

end

function vecintrot = reordint( vecintrot )

vecvecs = vecintrot - mean( vecintrot, 1) ;

theta = zeros( size( vecintrot, 1 ), 1 );

for i=2:size(theta,1)
    
    theta(i) = acos( dot( vecvecs(1,:), vecvecs(i,:) ) / ( norm(vecvecs(1,:)) * norm(vecvecs(i,:)) ) );
    
    n = cross( vecvecs(1,:), vecvecs(i,:) );
    
    if n(2)<0
        
        theta(i) = 2*pi - theta(i);
        
    end
    
end

[ ~, idx ] = sort(theta);

vecintrot = vecintrot( idx, : );

end