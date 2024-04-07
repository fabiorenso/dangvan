% MIT License
% 
% Copyright (c) 2023 fabiorenso
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% References:
% Dang Van, K. Multiaxial Fatigue Limit, a new approach. Preprint at (1984).
% Dang Van, K., Griveau, B. & Message, O. On a New Multiaxial Fatigue Limit Criterion Theory and Application. Biaxial and Multiaxial Fatigue 479–496 (1989).
% Dang Van, K., Cailletaud, G., Flavenot, J. F., Douaron, A. Le & Lieurade, H. P. Criterion for High Cycle Fatigue Failure Under Multuiaxial Loading. Biaxial and Multiaxial Fatigue 459–478 (1989).
% Lee, Y.-L., Barkey, M. E. & Kang, H.-T. Metal Fatigue Analysis Handbook: Practical Problem-Solving Techniques for Computer-Aided Engineering. (2012).

clear
clc
close all

%fatigue limit for reversed bending of the material
sigma_cr_f_inv=700.0; %[MPa]

%number of steps
t_steps=3;

%just initialize the stress to allocate the memory
stress.stress(1:6,1:t_steps)=0.0;
%S11 S22 S33 S12 S23 S31

%Create random stresses
% for n=1:t_steps
%     stress.stress(:,n)=rand(6,1)*100;
% end

%Or put your real stress
stress.stress(:,1)=[    0;    0;    0;    0;    0;    0];
stress.stress(:,2)=[  900;    0;    0;    0;  180;    0];
stress.stress(:,3)=[  900;    0;    0;    0; -180;    0];

%Dang Van parameters
a_dang=0.232;
b_dang=0.57735*sigma_cr_f_inv;
%--------------------------------------------------------------------------
% Real Dang Van parameters:									                                    
% a_dang = 3*((tau_cr_f_inv/sigma_cr_inv)-0.5)				                
% b_dang = tau_cr_inv [MPa]
%
% Approximated Dang Van parameters (considering Von Mises ratio between normal and shear stresses):
% a_dang: 3*(0.57735*sigma_cr_f_inv/sigma_cr_f_inv-0.5)=0.232
% b_dang: 0.57735*sigma_cr_f_inv [MPa]
%--------------------------------------------------------------------------

%just initialize all the variables to allocate the memory
stress.dev_res(1:6,1)=0;
stress.tau(1:t_steps)=0;
safety(1:t_steps)=0;
stress.idro(1:t_steps)=0;
r(1)=0;
l=1;
k=1;
rho(k)=r(l);
epsilon=1;
error=0.0001;
hardening_c=0.05;

%Compute hydrostatic and deviatoric stress for all timestep
%Start the matrix of deviatoric residual stress according to the paper
for j=1:t_steps
    stress.idro(j)=mean(stress.stress(1:3,j));
    stress.dev(:,j)=stress.stress(:,j)-stress.idro(j)*[1, 1, 1, 0, 0, 0]';
    stress.dev_res(:,1)=stress.dev_res(:,1)+stress.dev(:,j)/t_steps;
end

%Compute the residual deviatoric stress according to:
while epsilon>error
    for j=1:t_steps
        m=j+1;
        if j==t_steps
            m=1;
        end
        local_tensor=(stress.dev(:,m)+stress.dev_res(:,l))';
        stress.local_tensor=[local_tensor(1),local_tensor(4),local_tensor(5);
                             local_tensor(4),local_tensor(2),local_tensor(6);
                             local_tensor(5),local_tensor(6),local_tensor(3)];
        %Calcolo tensioni principali
        stress.princ=eig(stress.local_tensor);
        dist=(norm(stress.princ));
        if dist>r(l)
            r(l+1)=r(l)+hardening_c*(dist-r(l));
            stress.dev_res(:,l+1)=stress.dev_res(:,l)-(dist-r(l+1))/dist*(stress.dev(:,m)+stress.dev_res(:,l));
            l=l+1;
        end
    end
    k=k+1;
    rho(k)=r(l);
    epsilon=(rho(k)-rho(k-1))/rho(k-1);
end

%Compute the safety factor considering the residual stresses
for j=1:t_steps
    local_tensor=stress.dev(:,j)+stress.dev_res(:,l);
    stress.dev_micro=[local_tensor(1),local_tensor(4),local_tensor(5);
                         local_tensor(4),local_tensor(2),local_tensor(6);
                         local_tensor(5),local_tensor(6),local_tensor(3)];
    stress.princ_local=eig(stress.dev_micro);
    stress.tau(j)=(max(stress.princ_local)-min(stress.princ_local))/2;
    
    %Correction for negative hydrostatic pressure
    if stress.idro(j)<0
        safety(j)=b_dang/(stress.tau(j)-a_dang*stress.idro(j));
    else
        safety(j)=b_dang/(stress.tau(j)+a_dang*stress.idro(j));
    end
    %If safety factor is higher than 100 it is clipped to 100
    safety(j)=min(safety(j),100);
end

%The minimum safety factor along the path will be the safety factor
coeff=min(safety)
