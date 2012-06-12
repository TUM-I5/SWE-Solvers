function[size] = create_nc_testphi(path)

%250 meters should be enough (we dont want "extreme"-problems)
maxWaterHeight = 250;
maxHeightVariation = 50;
%300km/h = 5000 m/s should be enough
maxWaveSpeed = 5000;
%25km/h = 416.6666 m/s;
maxSpeedVariation = 416.666666;
dryTol = 0.001;
zeroTol = 0.0000001;

%delete old file
delete(path)

%ask for size
wallSize = input('Please enter the wall sample size: ');
fprintf('wall sample size: %d\n', wallSize);

randomSize = input('Please enter the random sample size: ');
fprintf('random sample size: %d\n', randomSize);

size = randomSize+wallSize;
fprintf('total sample size: %d\n', size);


%create variables
nccreate(path, 'hLow', ...
         'Dimensions' , {'sample', size}, ...
         'Format', 'classic' ... %netcdf4 seems to have some sort of limit: HDFError
         );

nccreate(path, 'hHigh', ...
         'Dimensions' , {'sample', size}...
        );

nccreate(path, 'huLow', ...
         'Dimensions' , {'sample', size}...
        );
    
nccreate(path, 'huHigh', ...
         'Dimensions' , {'sample', size}...
        );

nccreate(path, 'hStar', ...
         'Dimensions' , {'sample', size}...
        );

%oprn file (low level)
ncid = netcdf.open(path, 'NC_WRITE');

disp('filling the file with random values')
tic
%set random values
  %fill hLow with random values
  randHeights = rand(1, size)*maxWaterHeight+dryTol;
  varid = netcdf.inqVarID(ncid,'hLow');
  %counting starts at 0 for low level access
  %and at 1 for hight level acess
  netcdf.putVar(ncid, varid, 0, size, randHeights);
  
  varid = netcdf.inqVarID(ncid,'hHigh');
  %fill with the same value for walls
  netcdf.putVar(ncid, varid, 0, wallSize, randHeights(1:wallSize));
  
  %random values for the rest
  randHeights = randHeights(wallSize+1:size) + ...
    (2*(rand(1, randomSize)-0.5))*maxHeightVariation;
  randHeights = max(randHeights, dryTol+zeroTol);
  
  netcdf.putVar(ncid, varid, wallSize, randomSize, randHeights)
  
  %fill huLow with random values
  randSpeeds = transpose((2*(rand(1, size)-0.5))*maxWaveSpeed);
  netcdf.sync(ncid); %sync before reading
  hLowVec = ncread(path, 'hLow', 1, size);
  %hu = h * u
  huLowVec = randSpeeds .* hLowVec;
  varid = netcdf.inqVarID(ncid,'huLow');
  netcdf.putVar(ncid, varid, 0, size, huLowVec);
  
  varid = netcdf.inqVarID(ncid,'huHigh');
  %fill with the same values + opposite sign for walls
  netcdf.putVar(ncid, varid, 0, wallSize, -huLowVec(1:wallSize));
  
  %random values for the rest
  randSpeeds = randSpeeds(wallSize+1:size) + ...
    transpose((2*(rand(1, randomSize)-0.5))*maxSpeedVariation);
  netcdf.sync(ncid); %sync before reading
  hHighVec = ncread(path, 'hLow', wallSize, randomSize);
  %hu = h * u
  huHighVec = randSpeeds .* hHighVec;
  %random values for the rest
  netcdf.putVar(ncid, varid, wallSize, randomSize, huHighVec);
  
netcdf.close(ncid)
toc

ncwriteatt(path,'/','title', 'Computed solutions for the middle state which arises in the Shallow Water Equations for homogeneous Riemann Problems')
ncwriteatt(path,'/','author','Alexander Breuer')
ncwriteatt(path,'/','history', ...
  strcat(...
    'MATLAB R2011b, shallowwater.m (',...
    'maxWaterHeight=', num2str(maxWaterHeight), ' maxHeightVariation=', num2str(maxHeightVariation), ' maxWaveSpeed=', num2str(maxWaveSpeed), ' maxSpeedVariation=', num2str(maxSpeedVariation), ' dryTol=', num2str(dryTol), ' zeroTol=', num2str(zeroTol)...
    ,')'...
))
ncwriteatt(path,'/','institution','Technische Universität München, Department of Informatics, Chair of Scientific Computing')
ncwriteatt(path,'/','source','roots of an algebraic function')
ncwriteatt(path,'/','references','http://www5.in.tum.de')


ncwriteatt(path,'/','creation_date',datestr(now))
     
ncdisp(path)

%ncread(path, 'hLow', 1, size)
%ncread(path, 'hHigh', 1, size)
%ncread(path, 'huLow', 1, size)
%ncread(path, 'huHigh', 1, size)