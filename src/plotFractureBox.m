nodes_p =  obj.SimulationInput.ReservoirProperties.UnstructuredGridData.vertex.coord;
scatter3(nodes_p(:,1),nodes_p(:,2),nodes_p(:,3))
hold on
frac_nodes_p = fractures.geometry.coords;
frac_nodes_p = [frac_nodes_p; frac_nodes_p(1,:)];
%frac_nodes_p = [obj.geometry.coords;obj.geometry.coords(1,:)];
plot3(frac_nodes_p(:,1),frac_nodes_p(:,2),frac_nodes_p(:,3))