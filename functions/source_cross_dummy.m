function [J]=source_cross_dummy(c,parameters)
K=parameters.Model.K;
R = voxel_roi_map;
K =pinv(K);
U_map =R*K;
Nw =parameters.Dimensions.Nw;
for i=1:Nw
    J(:,:,i) = U_map*c(:,:,i)*U_map';
end
end
