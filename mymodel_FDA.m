function  [X] =  mymodel_FDA(HSI, MSI, R, FBm)

 [nr, nc,~]=size(MSI);
 L=size(HSI,3);
 sf = size(MSI,1)/size(HSI,1);
%  total_patches = 25;
total_patches = 100;
 rho1 = 0.0001;
 rho2 = 1e-4;
lambda = 1e-3;
 k = L;% PU:-5
 HSI_int=zeros(nr,nc,L);
 HSI_int(1:sf:end,1:sf:end,:)=HSI;

 %% Initiation
 FBmC  = conj(FBm);
 FBs  = repmat(FBm,[1 1 L]);   
 FBCs1=repmat(FBmC,[1 1 L]);
 HHH=ifft2((fft2(HSI_int).*FBCs1));
 HHH1=hyperConvert2D(HHH);
 
Zt=HSI_int;
HR_load1=imresize(HSI, sf,'bicubic');

ZE=hyperConvert2D(HR_load1);

MSI3=Unfold(MSI,size(MSI),3);

n_dr=nr/sf;
n_dc=nc/sf;

W=zeros([k,nr*nc]);

G1=zeros(size(ZE));
G2=zeros(size(W));
phi=zeros([L,k]);
Z_ba=zeros(size(W));

CCC=R'*MSI3+HHH1; 
C1=R'*R+(rho1/2)*speye(size(R,2)); 

knots = linspace(1, L, k); 
x = (1:L)';
for j = 1:k
    coeff = zeros(1, k); 
    coeff(j) = 1; 
    phi(:, j) = spline(knots, coeff, x); 
end
M1=rho1*(phi'*phi)+rho2*speye(size(phi,2));

%% ERS for suppixel
Y1=normalize(MSI3');
B=sparsepca(Y1);
Y1=Y1*B;
Y1=reshape(Y1,nr,nc);
[lables,~]=suppixel(Y1,total_patches);

%%  The ADMM-based solver
for i=1:50
ZE_last = ZE;


Z_ba = phi \ ZE; 

%Z
HR_HSI3=(rho1/2)*(phi*Z_ba)-(G1/2);

C3=CCC+HR_HSI3;
[ZE] = Sylvester(C1,FBm, sf,n_dr,n_dc,C3);  
Zt=hyperConvert3D(ZE,nr, nc );

%%Z_ba
M2=rho1*phi'*ZE+phi'*G1+rho2*W+G2;
Z_ba=M1 \ M2;

%%W
    for ll=0:1:total_patches-1
        index=find(lables==ll);
        ggg=Z_ba-G2/rho2;

        temp=ggg(:,index);
        C=lambda/rho2;
%      aa5= prox_nuclear (temp,C);
  
        aa3  =   repmat(mean( temp, 2 ),1,size(temp,2));
        aa4    =   temp - aa3;    
        aa5= prox_nuclear (aa4,C)+ aa3;
        W(:,index)=aa5;
    end
    
%%Gi
G1 = G1+rho1*(ZE-phi*Z_ba);
G2 = G2+rho2*(W-Z_ba);
  
   err(i) = norm(ZE_last(:) - ZE(:)) / norm(ZE_last(:));
    if  err(i)<1e-3
        break;
    end

end
X = Zt;
end