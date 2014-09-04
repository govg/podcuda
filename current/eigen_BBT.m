clc ;
clear all ;
format long e
load P_MATRIX_BBT.dat;
dim=size(P_MATRIX_BBT,1);
dim=dim^0.5;
bbt1=reshape(P_MATRIX_BBT,dim,dim);
[V2,D2]=eig(bbt1);
V2=reshape(V2,dim*dim,1);
fid =fopen('eigenvalue2.dat','w');
fid2 =fopen('eigenvector2.dat','w');
for i=1:dim%*dim%/2
   % for i=1:70%/2
        fprintf(fid,'%25.16f \n',D2(i,i));
        
   %last one f34.32
    %end
end
 
for i=1:dim*dim%1236544
  fprintf(fid2,'%25.16f \n',V2(i));
end
