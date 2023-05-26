for caty = 0.4:0.1:0.6
    %Establish pa   rameters for the simulation
    v   = caty;
    N   = 500;
    M   = uint16(N/v);
    c   = 1.0;
    s_c = sqrt(0.1);
    k   = 20.0;
    s_k = 0.0;
    m   = 1.0;
    au  = 1.0;
    rep = 20;
    
    %Establish  parameters for time
    tsim = 500;
    h    = 0.001;
    ite  = tsim/h;
    
    phi1  = zeros([rep,1],"double");
    psi1  = zeros([rep,1],"double");
    n_pob = zeros([N,rep],"double");
       
    for j=1:rep  
        %Declarate Necesary Matrix's
        ru   = ones([M,1],"double")*2.0;
        ru_v = zeros([M,1],"double");
        ni   = ones([N,1],"double")*10.0;
        ni_v = zeros([N,1],"double");
        mi   = ones([N,1],"double")*m;
        k1i  = zeros([N,1],"double"); k2i = zeros([N,1],"double"); k3i = zeros([N,1],"double"); k4i = zeros([N,1],"double");
        k1u  = zeros([M,1],"double"); k2u = zeros([M,1],"double"); k3u = zeros([M,1],"double"); k4u = zeros([M,1],"double");
        
        %rng('default') % For reproducibility
        ku   = normrnd(k,s_k,[M,1]);
        cu   = normrnd(c/N,s_c/sqrt(N),[N,M]);
        ci   = cu.';
        
        t    = 0.0;
        for i = 1:ite
            k1i = fi(ni,ru,mi,cu,t);
            k2i = fi(ni+h*k1i/2,ru,mi,cu,t+h/2);
            k3i = fi(ni+h*k2i/2,ru,mi,cu,t+h/2);
            k4i = fi(ni+h*k3i,ru,mi,cu,t+h);
        
            k1u = fu(ni,ru,ku,ci,au,t);
            k2u = fu(ni,ru+h*k1u/2,ku,ci,au,t+h/2);
            k3u = fu(ni,ru+h*k2u/2,ku,ci,au,t+h/2);
            k4u = fu(ni,ru+h*k3u,ku,ci,au,t+h);
        
            ni_v = ni + h*(k1i+2*k2i+2*k3i+k4i)/6;
            ru_v = ru + h*(k1u+2*k2u+2*k3u+k4u)/6;
        
            ni = ni_v;
            ru = ru_v;
        
            t = t + h;
        end
        
        posArrCount1=0;
        for i=1:N
            if ni(i,1) > 1/100000
                posArrCount1 = posArrCount1 + 1;
            end
        end
        
        posArrCount2=0;
        for i=1:M
            if ru(i,1) > 1/100000
                posArrCount2 = posArrCount2 + 1;
            end
        end
        
        phi1(j,1) = posArrCount1/N;
        psi1(j,1) = posArrCount2;
        n_pob(:,j)= ni;
        writematrix(phi1,'MacArthur Runge–Kutta phi_s m '+string(m)+';a '+string(au)+';c '+string(c)+';s_c^2 '+string(s_c^2)+';k '+string(k)+';s_k '+string(s_k)+';v '+string(v)+';N '+string(N)+'.csv')
        writematrix(psi1,'MacArthur Runge–Kutta psi_s m '+string(m)+';a '+string(au)+';c '+string(c)+';s_c^2 '+string(s_c^2)+';k '+string(k)+';s_k '+string(s_k)+';v '+string(v)+';N '+string(N)+'.csv')
        writematrix(n_pob,'MacArthur Runge–Kutta n_pob m '+string(m)+';a '+string(au)+';c '+string(c)+';s_c^2 '+string(s_c^2)+';k '+string(k)+';s_k '+string(s_k)+';v '+string(v)+';N '+string(N)+'.csv')
        fprintf('Iteración: '+string(j)+'\n')
     end
end

function fit = fi(ni, ru, mi, cu, ~)
    fit = cu*ru-mi;
    fit = fit.*ni;
end

function fut = fu(ni, ru, ku, ci, au, ~)
    fut = au*(ku-ru)-ci*ni;
    fut = fut.*ru;
end
