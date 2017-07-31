%subplot(1,3,1)
Y1 = mdsify(mean(sim_struct_indv,3),'rtoz');
%text(Y1(:,1)+0.01,Y1(:,2),conditions)


%subplot(1,3,2)
Y2 = mdsify(models(1).rsa,'rtoz');
%text(Y2(:,1)+0.01,Y2(:,2),conditions)

%subplot(1,3,3)
[D,Z,T] = procrustes(Y1,Y2,'Scaling',true);
plot(Y1(:,1),Y1(:,2),'bo',Z(:,1),Z(:,2),'rd');
text(Y1(:,1)+0.01,Y1(:,2),conditions)
text(Z(:,1)+0.01,Z(:,2),conditions)
