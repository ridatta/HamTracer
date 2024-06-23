classdef optics
    properties
    end
    methods
        function obj = optics()            
        end
        
        function out = get_ref_ang(~,thi, n1, n2) % snells law
            % thi = incidence angle [rad]
            out = asin( n1./n2 .* sin(thi)); 
        end
        
        % Fresnel coefficients
        function out = get_r(obj,n1,n2,thi,pol) 
            tht = obj.get_ref_ang(thi, n1, n2);
            if strcmpi(pol,'s')
                out = (n1 .* cos(thi) - n2 .* cos(tht)) ./ (n1 .* cos(thi) + n2 .* cos(tht));
            elseif strcmpi(pol,'p')
                out = (n2 .* cos(thi) - n1 .* cos(tht)) ./ (n2 .* cos(thi) + n1 .* cos(tht));
            end
        end
     
         function out = get_t(obj,n1,n2,thi,pol)
             tht = obj.get_ref_ang(thi, n1, n2);
            if strcmpi(pol,'s')
                out = 2 * n1 .* cos(thi) ./ (n1 .* cos(thi) + n2 .* cos(tht));
            elseif strcmpi(pol,'p')
                out = 2 * n1 .* cos(thi) ./ (n2 .* cos(thi) + n1 .* cos(tht));
            end
         end
         
         function out = getReflectance(obj,n1,n2,thi,pol)
             r = obj.get_r(n1,n2,thi,pol);
             out = r.^2; 
         end
         
          function out = getTransmittance(obj,n1,n2,thi,pol)
             tht = obj.get_ref_ang(thi, n1, n2);
             t = obj.get_t(n1,n2,thi,pol);
             out = (n2 .* cos(tht)) ./ (n1 .* cos(thi)) .* t.^2;
          end
         
          function out = getBrewsterAng(~,n1,n2) % Brewster angle, R_p = 0
              out = atan(n2./n1);
          end
          
          function out = getCritAng(~,n1,n2) % critical angle, R_p = 0
              out = asin(n2./n1);
          end
          
          function out = getThinFilmReflectance(obj,n1,n2,n3,d,thi,pol,lam)
              % [n1, n2, n3] = refractive idx. of layers
              % d = thickness [m]
              % thi = incidence nagle in medium 1, rad
              % pol = 's' or 'p' polarization
              % lam = [m] wavelength
              
                k0 = 2 * pi / lam;  % [m] 
              
                r12 = -1 * obj.get_r(n1,n2,thi,pol); 
                th2 =  obj.get_ref_ang(thi, n1, n2); 

                r23 = -1* obj.get_r(n2,n3,th2,pol);

                phi = 1 * n2 * k0 * d * cos(th2); 

                r = (r12 + r23 * exp(1i * 2 * phi)) ./ (1 + r12 * r23 * exp(1i * 2 * phi)); 
                out = r .* conj(r);
          end
              
          
          % farby-Perot
          function out = getDel(~,k,n2,L,tht)
              out = k * n2 * L / cos(tht);
          end
          function out = getFP_r(obj,n1,n2,thi,pol,k,L)
              tht = obj.get_ref_ang(thi, n1, n2);
              del = obj.getDel(k,n2,L,tht);
              r12 = obj.get_r(n1,n2,thi,pol);
              out = r12 .* (1 - exp(2*1i*del)) ./ (1 - r12.^2 * exp(2*1i*del)); 
          end
            function out = getFP_t(obj,n1,n2,thi,pol,k,L)
              tht = obj.get_ref_ang(thi, n1, n2);
              del = obj.getDel(k,n2,L,tht);
              r12 = obj.get_r(n1,n2,thi,pol);
              out = (1 - r12.^2) * exp(1i*del) ./ (1 - r12.^2 * exp(2*1i*del));  
            end
          
            function out = getFPTransmission(~,R,del)
                out = (1 - R.^2) ./ ( (1-R).^2 + 4 .* R .* sin(del).^2);
            end
            
            function out = getFPTReflection(~,R,del)
                out = 4 * R .* sin(del).^2 ./ ( (1-R).^2 + 4 * R .* sin(del).^2);
            end
          
          
         
    end
end