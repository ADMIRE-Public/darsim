function avgDist_rock = Obtain_Average_Distance_Plane_From_Tetrahedron( intersectCoord, TetrahedralGridData, DistPontInTetra, Im )
%

TetraCoord = TetrahedralGridData.Vertex.Coordinates( TetrahedralGridData.Cell.Vertices(Im,:) , : );

MeioInt = mean( intersectCoord, 1 );

PontosTetra = DistPontInTetra*TetraCoord;

VecTetra = PontosTetra - repmat( MeioInt, [size(PontosTetra,1) 1] );

if size(intersectCoord,1)>1

    NormInt = cross( intersectCoord(1,:) - MeioInt, intersectCoord(end-1,:) - MeioInt ); NormInt = repmat( NormInt/norm(NormInt), [size(VecTetra,1) 1] );

    dist = abs( sum( VecTetra.*NormInt, 2) );

    % Determina  a distância média entre a fratura e o tetraedro ----------
    avgDist_rock = sum( dist ) / size( dist, 1 );
    %----------------------------------------------------------------------

else
    
    % Determina  a distância média entre a fratura e o tetraedro ----------
    avgDist_rock = 1;
    %----------------------------------------------------------------------
    
end

end