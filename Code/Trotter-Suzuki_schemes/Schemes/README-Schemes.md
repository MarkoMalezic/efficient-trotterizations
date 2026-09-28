# Trotter-Suzuki decompostion schemes data

In this directory, we collect the scheme data from our minimization attempts, gathered inside `Trotter-Suzuki/`.
The schemes are first structured by their order and then the number of cycles they require.

## All schemes

In folder `Trotter-Suzuki/All_schemes` we archive schemes at orders n = 2, 4, 6, and 8 for number of cycles q < 31.
We store three different precision types (double, long double and quad).
At double precision we have ~ 15 significant digits, long double expands this to ~ 18 significant digits, and quad doubles it to ~ 34 significant digits.

For double and long double, each cycle has its own `cycles_q#` folder, which firstly includes the raw data file `schemes_q#.out`.
The scheme parameters inside this file were minimized twice, once to find the minimum region and again to pinpoint the exact values to numerical precision.
An important note here is that these scheme parameters are multiplied by 1/2 to those found in Section 2.2 of the paper "Efficient Trotter-Suzuki Schemes for Long-Time Quantum Dynamics" (see Eq. (10)).

Secondly, these scheme parameters and their errors/efficiencies were then transformed into Numpy arrays in order to be analyzed with Python.
Additionally to the scheme parameters in the notation of our implementation, we transform them into a "standard" notation, which is the one used by Omelyan et al. (see Eq. (8)).
We also transform these into the general notation with the c and d scheme parameters as defined by Eq. (9).
From these we calculate the distance from the origin (uniformity), and also save these values.

Finally, we collect all of these arrays in an HDF5 format in case they need to be read differently. However, loss of precision is possible due to HDF5 not supporting long_double/float128 in all formats.

Quad precision only provides the raw `schemes_q#.out` file for each cycle (placed directly under `All_schemes/quad/order_n#/`, without a `cycles_q#` subfolder), since the notation transforms, distance-from-origin calculation and HDF5 export are not implemented for quad precision.

## Recommended schemes

In the second folder `Trotter-Suzuki/Recommended_schemes`, we select specific schemes from our archive, which are believed to perform well in practice due to the theoretical efficiency and uniformity arguments.
We do this at each no. cycles q, and in each `cycles_q#` folder you can find `schemes_q#.out`, which contains the recommended scheme with a brief description of their properties.
As before we transform them into other notations, compute the properties, and package them in Numpy arrays and HDF5 format for double and long double precision, each stored in its own `double`/`long_double` subfolder.
