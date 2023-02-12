# LABOGRID

An automatic ligand-based gridbox size calculation tool for molecular docking.

```
Usage
╰─○ labogrid.py [-l] <filename>

Argument
╰─○ Command description:
        -l  ligand filename (supported: pdb, pdbqt, sdf, mol2)
        -e  experimental ligand filename (supported: pdb, sdf, mol2)
        -h  help
        -a  about
╰─○ Optional parameters:
        -s  scale factor (default is 2)

Return
╰─○ Gridbox Center: X  {value} Y  {value} Z  {value}
    Gridbox Size  : W  {value} H  {value} D  {value}
```

## Method
1. `min()` and `max()` were used to determine the minimum and maximum X, Y, Z atomic coordinate of a ligand.
2. Using the values in 1, `statistics.mean()` was used to determine the X, Y, Z coordinate of a gridbox center.
3. Using the values in 1, the `abs(value)` of subtraction between the minum and maximum X, Y, Z atomic coordinate of a ligand was used to determine the size of gridbox in terms of width, length and depth. 
4. **Scale factor of 2** (default) was used to adjust the gridbox size based on values in 3.

> **NOTE:** All values are rounded to 3 decimals at the end of calculation. Further testing will be performed to determine the optimal scale factor. 

## Bug
If you encounter any bugs, please report the issue to https://github.com/RyanZR/labogrid/issues.

## License
The script is licensed under MIT, see the [LICENSE](https://github.com/RyanZR/labogrid/blob/main/LICENSE) file for details.
