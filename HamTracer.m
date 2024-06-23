classdef HamTracer
    properties
        lam = 500e-9; % [m], free space light wavelength
        c = 3 * 1e8; % m/s
        w; % rad/s
        n_grid = 101; % grid points
        xx; % x coordinates
        zz; % z coodinates
        dx; % [m] grid size
        N; % refractive idx
        dxN; 
        dzN;
        results = struct();
    end
    methods
        function obj = HamTracer(lam,n_grid,xlimits,zlimits)
            % n_grid = grid size in x-dimension
            % xlimits = [xmin, xmax] in [m]
            % zlimits = [zmin, zmax] in [m]
            obj.lam = lam;
            obj.w = 2 * pi / lam * obj.c; % rad / s
            obj.n_grid = n_grid;
            x =  1 * linspace(xlimits(1),xlimits(2),n_grid); % Grid, [m]
            obj.dx = max(diff(x));
            z =  zlimits(1):obj.dx:zlimits(2); % [m]
            [obj.xx,obj.zz] = meshgrid(x,z); % create grid
        end
        function out = getOwaveN(obj,ne)
            addpath('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Codes/Other/PlasmaFormulary/')
            addpath('C:\Users\rdatta\Dropbox (MIT)\PUFFIN\Codes\Other\PlasmaFormulary\')
            F = FundamentalPlasma();
            wp = F.getElecPlasmaFreq(ne); % rad/s
            out = sqrt(1 - wp.^2/obj.w.^2); % [-]
        end
        function out = getXwaveN(obj,ne,B)
            addpath('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Codes/Other/PlasmaFormulary/')
            addpath('C:\Users\rdatta\Dropbox (MIT)\PUFFIN\Codes\Other\PlasmaFormulary\')
            F = FundamentalPlasma();
            wp = F.getElecPlasmaFreq(ne); % rad/s
            wh = F.getElecPlasmaFreq(ne) + F.getElectronGyroFreq(SI2cgs(B,'Magnetic Field')); % upper hybrid freq. rad/s
            out = sqrt(1 - wp.^2/obj.w.^2 .* (obj.w^2 - wp.^2) ./ (obj.w^2 - wh.^2));
        end
        function obj = setN(obj,N)
            % Sets the n_grid refractive index 
            obj.N = N;
           [obj.dxN,obj.dzN] = gradient(N,obj.dx,obj.dx); % x- and z-gradient of N
        end
        function obj = solve(obj,xi,zi,thi,tend,tstep)
            % xi = initial x-position of rays in [m]
            % thi = initial angle of each ray
            % tspan = stopping time
            % tstep = max. step size
            T = {}; xout = {}; zout = {}; kxout = {}; kzout = {};
            for ii = 1:numel(xi)
                    x0 = xi(ii);
                    z0 = zi(ii);
                    k0 = interp2(obj.xx,obj.zz,obj.N,x0,z0) * obj.w / obj.c; % initial wavenumber
                    th0 = thi(ii); % Initial angle
                    kx0 = k0 * sin(th0);
                    kz0 = k0 * cos(th0);

                    y0 = [x0,z0,kx0,kz0]; % IC
                    [t{ii},yout] = ode45(@(t,y) obj.ode(t,y),[0,tend],y0,...
                        odeset('MaxStep',tstep,'Events',@obj.myEvent));
                    xout{ii} = yout(:,1);
                    zout{ii} =  yout(:,2);
                    kxout{ii} =  yout(:,3);
                    kzout{ii} = yout(:,4);
            end
            obj.results.t = t;
            obj.results.xout = xout;
            obj.results.zout = zout;
            obj.results.kxout = kxout;
            obj.results.kzout = kzout; 
            
        end
        function out = ode(obj,t,y) % ODE to solve
          out = zeros(4,1);

          x = y(1); 
          z = y(2);
          kx = y(3);
          kz = y(4);
            
          n = interp2(obj.xx,obj.zz,obj.N,x,z); % n at current (x,z)
          dxn = interp2(obj.xx,obj.zz,obj.dxN,x,z); % n at current (x,z)
          dzn = interp2(obj.xx,obj.zz,obj.dzN,x,z); % n at current (x,z)

          if isnan(n)
              n = interp2(obj.xx,obj.zz,obj.N,max(obj.xx(:)),max(obj.zz(:)));
          end
          if isnan(dxn)
              dxn = 0;
          end
          if isnan(dzn)
              dxn = 0;
          end

          dHdkx = 2 * kx;
          dHdkz = 2 * kz;
          dxH = -obj.w^2/obj.c^2 * 2 * n * dxn;
          dzH = -obj.w^2/obj.c^2 * 2 * n * dzn;

          dxdt = dHdkx;
          dzdt = dHdkz;
          dkxdt = -1 * dxH;
          dkzdt = -1 * dzH;

          out(1) = dxdt;
          out(2) = dzdt;
          out(3) = dkxdt;
          out(4) = dkzdt;
        end
        function [value, isterminal, direction] = myEvent(obj,t, y)
                value      = (y(2) >=  0.98*max(obj.zz(:)));
                isterminal = 1;   % Stop the integration
                direction  = 0;
        end
        function showN(obj)
            surf(obj.zz * 1e3,obj.xx * 1e3,obj.N); colormap('gray'); colorbar();
            xlabel('z [mm]');
            ylabel('x [mm]')
            view(2);
            title('Refractive Index');
            shading interp
            formatPlots(600);legend('off');
            axis equal
            ylim(1e3*[min(obj.xx(:)),max(obj.xx(:))])
            xlim(1e3*[min(obj.zz(:)),max(obj.zz(:))])
        end
        
        function showResults(obj,varargin)
            if nargin > 1
                clr = varargin{1};
            else
                clr = 'r';
            end
            for ii = 1:numel(obj.results.t)
                plot(obj.results.zout{ii}*1e3,...
                    obj.results.xout{ii}*1e3,...
                    '-','Color',clr,'Linewidth',2); hold on;
            end

            xlabel('z [mm]')
            ylabel('x [mm]')
            grid on;
            formatPlots(600); legend('off')
            axis equal
            ylim(1e3*[min(obj.xx(:)),max(obj.xx(:))])
            xlim(1e3*[min(obj.zz(:)),max(obj.zz(:))])
            
        end

    end
end