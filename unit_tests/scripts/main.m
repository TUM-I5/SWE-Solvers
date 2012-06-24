% Author: Alexander Breuer, breuera AT in.tum.de
%
% This function creates a .nc files, fills it with random values and
% computes the hStar values (see create_nc_testphi)
%
function[] = main()

path = input('Please enter the path (including the filename),\nwhere the file should be written to (needs to be quoted ''/home/[...]/tesphi.nc''): ');
fprintf('file path: %s\n', path);

size = create_nc_testphi(path)

hLow = ncread(path, 'hLow', 1, size);
hHigh = ncread(path, 'hHigh', 1, size);
huLow = ncread(path, 'huLow', 1, size);
huHigh = ncread(path, 'huHigh', 1, size);
hStar = zeros(size, 1);

disp('computing middle heights');
tic
for j=1:size
  hStar(j) = calculate_hstar(hLow(j), hHigh(j), huLow(j), huHigh(j));
  if mod(j ,(size/10)) == 0
      fprintf('%i %% done.\n', (j*100/size));
  end
end
toc

disp('writing hStar to the file')

%open file (low level)
ncid = netcdf.open(path, 'NC_WRITE');
varid = netcdf.inqVarID(ncid,'hStar');
%write the values
netcdf.putVar(ncid, varid, 0, size, hStar);
netcdf.close(ncid)

disp('script finished');