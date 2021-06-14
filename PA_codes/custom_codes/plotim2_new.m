function out = plotim2_new(mesh,val)
if mesh.dimension == 3 && ...
        (strcmp(mesh.type,'stnd_bem') || strcmp(mesh.type,'spec_bem') || strcmp(mesh.type,'fluor_bem'))
    h = trisurf(mesh.elements,...
	    mesh.nodes(:,1),...
	    mesh.nodes(:,2),...
	    mesh.nodes(:,3),...
	    val,'FaceAlpha',0.5);
else
    h = trisurf(mesh.elements,...
	    mesh.nodes(:,1),...
	    mesh.nodes(:,2),...
	    mesh.nodes(:,3),...
	    val);
       
%     out = zeros(351,501);
    out = zeros(601,601);
    for i=1:length(h.Vertices(:,1))
%         out(floor((h.Vertices(i,2)*5+1)),floor(((h.Vertices(i,1)+50)*5+1)))=h.FaceVertexCData(i);%Change the 50 if mesh changes
        %out(round((h.Vertices(i,2)*5+1)),round(((h.Vertices(i,1)+50)*5+1)))=h.FaceVertexCData(i);%Change the 50 if mesh changes
        out(round((h.Vertices(i,2)*10+1)),round(((h.Vertices(i,1)+30)*10+1)))=h.FaceVertexCData(i);
    end
    out=flip(out);    
end
shading interp;
view(2);
axis equal; 
axis tight
axis on;
colormap hot;