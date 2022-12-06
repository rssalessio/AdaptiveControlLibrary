%Compute L1 norm of a MIMO system
% Slow implementation
%[L1norm] = L1norm(SYS);
%[L1norm] = L1norm(SYS,tollerance);
function [L1] = L1norm_mimo(SYS,toll)
%% Setting up
if (nargin == 1)
    toll = 10^(-2);%default tollerance for impulse response
end
%% Verifying the model
[zdata,pdata,kdata] = zpkdata(zpk(SYS));
MosdelDim  = size(kdata);
%checking gain case
gf = 1;
for i=1:MosdelDim(1)
    for j=1:MosdelDim(2)
        if (size(pdata{i,j},1)~=0)
            gf = 0;
        end
    end
end
if (gf==1)
    L1i = zeros(MosdelDim(2),1);
    for i=1:MosdelDim(1)
        for j=1:MosdelDim(2)
            L1i(i) = L1i(i) + abs(kdata(i,j));
        end
    end
    L1 = max(L1i);
    return
end
%Checking stability
maxabsp = 0;
for i=1:MosdelDim(1)
    for j=1:MosdelDim(2)
        maxabspcand = max(abs(pdata{i,j}));
        if (maxabsp < maxabspcand)
            maxabsp = maxabspcand;
        end
        if (max(pdata{i,j})>=0)
            L1 = inf;
            warning('EK:L1norm:sytem_unstable',...
                'The system is (internally) unstable');
            return
        end
    end
end
%checking properness
for i=1:MosdelDim(1)
    for j=1:MosdelDim(2)
        if (size(zdata{i,j},1)>size(pdata{i,j},1))
            warning('EK:L1norm:sytem_improper','The system is improper');
            L1 = inf;
            return
        end
    end
end
%% Some needed preparing computations
SYSss = ss(SYS);
%% Nested functions
%defining internal values
L1c = zeros(size(SYSss.c,1),1);
Stopf = false;
xf = zeros(size(SYSss.b,1)+size(SYSss.c,1),1);
%Defining a nested function for solver (dx/dt)
    function dx = dxdt(t,x)
        %extracting the system's state
        xr = x(1:size(SYSss.c,2));
        %Computing the state derivatives
        dxr = SYSss.a*xr;
        %Computing the impulse integral
        dxi = abs(SYSss.c*xr);
        %Combining the derivatives
        dx = [dxr;dxi];
    end
%Defining a nested output function
    function status = OutFcn(t,x,flag)
        if ~(strcmp(flag,'init')||(strcmp(flag,'done')))
           %extracting the system's state
            xe = x(:,size(x,2));
            xr = xe(1:size(SYSss.c,2));
            %Computing the norm of state
            xn = norm(xr);
            %Making decision
            if (xn<=toll)
                L1c = xe((size(SYSss.c,2)+1):(size(xe,1)));
                Stopf = true;
                status = 1;
            else
                xf = xe;
                status = 0;
            end
        else
            status = 0;
        end
    end
%% Computing the L1 norm matrix
%Setting Solver properties
options = odeset('MaxStep',(1/(2*maxabsp)),'AbsTol',toll,'OutputFcn',@OutFcn);
%options = odeset('AbsTol',toll,'OutputFcn',@OutFcn);
%Computing span of t (not niether step nor final time)
dt = min(abs(pdata{1,1}));
for i=1:MosdelDim(1)
    for j=1:MosdelDim(2)
        dt = min(dt,min(abs(pdata{i,j})));
    end
end
dt = 1/dt;
%INTEGRATION CICLE
L1v = zeros(size(SYSss.c,1),1);
for j=1:MosdelDim(2)
    %Setting initial conditions
    x0 = [SYSss.b(:,j);zeros(size(SYSss.c,1),1)];
    while ~Stopf
        ode45(@dxdt,[0 dt],x0,options);
        x0 = xf;
    end
    L1v = L1v + L1c;
end
%Computing the L1 norm     
L1 = max(L1v);
end