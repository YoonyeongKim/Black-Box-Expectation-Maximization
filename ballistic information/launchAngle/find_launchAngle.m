clear all; clc;

site1BM3 =  importdata('site1BM1.mat') ;
[X,Y] = meshgrid(site1BM3.AzDeg,site1BM3.ElDeg);
Z1 = site1BM3.Latitude;
Z2 = site1BM3.Longitude;


[ans1, dummy] = contour(X,Y,Z1-37.5,[0 0]);
hold on
[ans2, dummy] = contour(X,Y,Z2-126.5,[0 0]);


for i = 1:size(ans1,2)
    for j = 1:size(ans2,2)
        dis(i,j) = norm(ans1(:,i)-ans2(:,j));
        I(i,j) = i; J(i,j) = j;
    end
end

dis_vec = dis(:); I_vec = I(:); J_vec = J(:);
[dummy, ind] = sort(dis_vec);

for ii = 1:10
    i = I_vec(ind(ii));
    j = J_vec(ind(ii));
    fa(:,ii) = ans1(:,i);
    text(fa(1,ii),fa(2,ii),num2str(ii))
end

