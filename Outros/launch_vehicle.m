function varargout=launch_vehicle(varargin)
%LAUNCH_VEHICLE Outputs appropriate launch vehicle for a mission
%
% [LV_DATABASE] = LAUNCH_VEHICLE
% Outputs a structure containing information on each launch vehicle.
%
% [LV_INDEX,APOGEE,PERIGEE] = LAUNCH_VEHICLE(RADIUS,MASS)
% Outputs the index into the LV_DATABASE of the launch vehicle that
% can most cost-effectively deliver a satellite of MASS kilograms
% into a circular parking orbit of RADIUS kilometers. If no launch
% vehicle can meet the demands, LV_INDEX is set to zero. APOGEE and
% PERIGEE are the actual values achievable for the specified MASS.
% Assumptions
% 1. Best launch vehicle is determined as a function of cost
% 2. Launch vehicle only needs to get payload to parking orbit
% 3. Reliability (i.e. 98% successful Space Shuttle) is not
% incorporated into the analysis
% 4. Lead time is not factored into the module
% 5. Launch occurs from Cape Canaveral (28.5 inclination)
% Input launch vehicle characteristics, including Excel performance equations
% 1999 AIAA International Reference Guide to Space Launch Systems
% Inflation factor for 1999 launch vehicle estimates (source: SMAD, Chapter 20)
% min and max are range of validity for performance equation
% 0 perigee means circular orbit
i=0;
i=i+1;lvdata(i) = struct('name','Atlas IIA',...
 'company','Lockheed Martin',...
 'cost',91,...
 'x2',NaN,'x',NaN,'intercept',NaN,... 
 'c_log',-847.55,'c_int',12041,...
 'perigee',185,...
 'min',0,'max',36000);

i=i+1;lvdata(i) = struct('name','Atlas IIAS',...
 'company','Lockheed Martin',...
 'cost',112,...
 'x2',NaN,'x',NaN,'intercept',NaN,...
 'c_log',-968.15,'c_int',13978,...
 'perigee',185,...
 'min',0,'max',36000);

i=i+1;lvdata(i) = struct('name','Atlas IIIA',...
 'company','Lockheed Martin',...
 'cost',112,...
 'x2',NaN,'x',NaN,'intercept',NaN,...
 'c_log',-969.02,'c_int',14319,...
 'perigee',185,...
 'min',0,'max',36000);

i=i+1;lvdata(i) = struct('name','Atlas IIIB',...
 'company','Lockheed Martin',...
 'cost',112,...
 'x2',NaN,'x',NaN,'intercept',NaN,...
 'c_log',-1248.6,'c_int',17808,...
 'perigee',185,...
 'min',0,'max',36000);
i=i+1;lvdata(i) = struct('name','Delta II 7320',...
 'company','Boeing',...
 'cost',59,... 
 'x2',0,'x',-1.6,'intercept',4800,...
 'c_log',NaN,'c_int',NaN,...
 'perigee',0,...
 'min',0,'max',2000);
i=i+1;lvdata(i) = struct('name','Delta II 7420',...
 'company','Boeing',...
 'cost',59,...
 'x2',0,'x',-1.4815,'intercept',4888.9,...
 'c_log',NaN,'c_int',NaN,...
 'perigee',0,...
 'min',0,'max',2000);

i=i+1;lvdata(i) = struct('name','Delta II 7920',...
 'company','Boeing',...
 'cost',64,...
 'x2',0,'x',-1.0526,'intercept',5789.5,...
 'c_log',NaN,'c_int',NaN,...
 'perigee',0,...
 'min',0,'max',2000);

i=i+1;lvdata(i) = struct('name','Delta III',...
 'company','Boeing',...
 'cost',96,...
 'x2',NaN,'x',NaN,'intercept',NaN,...
 'c_log',-848.91,'c_int',12723,...
 'perigee',0,...
 'min',0,'max',36000);

i=i+1;lvdata(i) = struct('name','Pegasus XL',...
 'company','Orbital Sciences',...
 'cost',16,... 
 'x2',0,'x',-4,'intercept',2000,...
 'c_log',NaN,'c_int',NaN,...
 'perigee',0,...
 'min',0,'max',2000);
i=i+1;lvdata(i) = struct('name','Space Shuttle',...
 'company','NASA',...
 'cost',320,...
 'x2',-.025,'x',-6,'intercept',25300,...
 'c_log',NaN,'c_int',NaN,...
 'perigee',0,...
 'min',0,'max',600);

i=i+1;lvdata(i) = struct('name','Taurus 2110',...
 'company','Orbital Sciences',...
 'cost',21,...
 'x2',0,'x',-2.1622,'intercept',3243.2,...
 'c_log',NaN,'c_int',NaN,...
 'perigee',0,...
 'min',0,'max',2000);

% determine which output type given the arguments
if (nargin==0 & nargout<=1)
 % return the LV databse
 varargout{1} = lvdata;
 return
elseif (nargin==2 & nargout==3)
 % eventually return the LV index
 A_po = varargin{1}-6378.136; % change radius to altitude
 mass = varargin{2};
else
 % something is wrong
 error('Invalid number of input or output arguments');
end
% Get matrices
cost = cat(2,lvdata.cost);
x2 = cat(2,lvdata.x2);
x = cat(2,lvdata.x);
c_log = cat(2,lvdata.c_log)
c_int = cat(2,lvdata.c_int);
intercept = cat(2,lvdata.intercept);
perigee = cat(2,lvdata.perigee);
amin = cat(2,lvdata.min);
amax = cat(2,lvdata.max);
% Eliminate launch vehicles unable to launch to desired apogee radius
bound_indices = find(A_po<=amax & A_po>=amin);
% Compute maximum payload capability for a given altitude
pcap_a = x2.*A_po.^2 + x.*A_po + intercept;
pcap_b = c_log.*log(A_po)+c_int;
for i=1:length(lvdata)
 if ~isnan(pcap_a(i))
 pcap(i) = pcap_a(i);
 elseif ~isnan(pcap_b(i))
 pcap(i) = pcap_b(i);
 else
 pcap(i) = 0;
 end
end
% Eliminate launch vehicles unable to meet mass requirement
cap_indices = find(pcap>=mass); 
% intersection of the sets
indices = intersect(cap_indices, bound_indices);
% find the minimum cost vehicle
[cost, index] = min(cost(indices));
index = indices(index);
if isempty(index)
 varargout{1} = NaN;
 varargout{2} = 0;
 varargout{3} = 0;
 return
end
if x2(index)==0
 apoapsis = (mass-intercept(index))/x(index);
elseif isnan(x2(index))
 apoapsis = min([exp((mass-c_int(index))/c_log(index)) amax(index)]);
else
 apoapsis = (-x2(index)-sqrt(x(index)^2-4*x2(index)*...
 (intercept(index)-mass)))/(2*x2(index));
end
if perigee==0
 periapsis=apoapsis;
else
 periapsis=perigee(index);
end
varargout{1} = index;
varargout{2} = apoapsis;
varargout{3} = periapsis; 
