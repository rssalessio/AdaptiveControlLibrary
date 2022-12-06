%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: l1norm(sys, tol, maxiter)          
% Author: Alessio Russo
% E-mail: russo.alessio890@gmail.com
% Date: 30/08/2016      
% Description: calculate L1 norm of state space system.
% Please use minreal to remove the zero/pole cancellation with an
% appropriate tolerance value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = l1norm(sys, tol, maxiter)
    if (isproper(sys) == 0)
        disp('The system is not proper');
        L=inf;
    end
    alpha=max(real(eig(sys)));
    if nargin<3, maxiter=20; end% maximum iterations
    if nargin<2, tol=-alpha*.001; end% tolerance
    
    sys=ss(sys);
    system_size = size(sys);
    L = zeros(system_size(1),1);
    
    for i = 1:system_size(1)
       for j = 1:system_size(2)
           red_sys = ss(sys.A, sys.B(:,j), sys.C(i,:), sys.D(i,j));
           %check if system is strictly proper. If the system is not
           %strictly proper we cannot use l1norm_siso but l1norm_mimo
           nz = length(zero(red_sys));
           np = length(pole(red_sys));
           if (nz==np)
               L(i)=L(i)+l1norm_mimo(red_sys, tol);
           else
               L(i)=L(i)+l1norm_siso(red_sys, tol, maxiter);
           end
       end
    end
    L = max(L);
end