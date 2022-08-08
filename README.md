# GridSphere

Julia implement of [GridSphere](
https://jp.mathworks.com/matlabcentral/fileexchange/28842-grid-sphere).

## Install
```julia
using Pkg
Pkg.add("GridSphere")
```

## Usgae
```julia
using GridSphere
k = 5
latlon = gridsphere(2 + 10 * 4 ^ k)
lats = latlon[:,1]
lons = latlon[:,2]
```

## Referece

[1] Kurt von Laven (2022). Grid Sphere (https://www.mathworks.com/matlabcentral/fileexchange/28842-grid-sphere), MATLAB Central File Exchange. Retrieved August 8, 2022.