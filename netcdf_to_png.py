import argparse
import netCDF4 as nc
import cv2
import numpy as np



parser = argparse.ArgumentParser(description="Convert NetCDF to png")
parser.add_argument('--input_file', help="Input Grib file", action='store', type=str)
parser.add_argument('--out_w', default='0', help="Output data width (0 = same as input)", action='store', type=int)
parser.add_argument('--out_h', default='0', help="Output data height (0 = same as input)", action='store', type=int)    
parser.add_argument('--output_file', default='output.png', help='Output file')
args = parser.parse_args()

print(f"Arguments: {args}")

ncdf_name = args.input_file
output_w = args.out_w
output_h = args.out_h
output_name = args.output_file

if (ncdf_name is None or ncdf_name == ""):
    print(f"Invalid arguments - missing file name")
    quit()

print("==========")
print(f"Loading NetCDF: {ncdf_name}")
        
    
ds = nc.Dataset(ncdf_name)  
pmsl = ds['msl'][:].squeeze()     

print("==========")
print("NetCDF Info")
print(f"Min values: {pmsl.min()}, Max value: {pmsl.max()}")


#convert values to range 0 - 255 (single byte)
values = np.interp(pmsl, (pmsl.min(), pmsl.max()), (0, 255))

#flip data vertically
#values = np.flipud(values)

#Auto-compute width or height if only one size is set
if (output_w == 0 and output_h != 0):    
    ar = values.shape[1] / values.shape[0]
    output_w = int(output_h * ar)

if (output_w != 0 and output_h == 0):    
    ar = values.shape[0] / values.shape[1]
    output_h = int(output_w * ar)

#resize if needed
if (output_w != 0 and output_h != 0):
    values = cv2.resize(values, (output_w, output_h))  # , interpolation=cv2.INTER_CUBIC)
else:
    output_w = int(values.shape[1])
    output_h = int(values.shape[0])

print("==========")
print("Output info")
print(f"Resolution: {output_w} x {output_h}")
print(f"File: {output_name}")

cv2.imwrite(output_name, values) 