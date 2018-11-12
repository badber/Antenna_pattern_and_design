Compute radiation patterns for different types of **broadside** and **end-fire** antenna-arrays.

(see [here](http://www.idc-online.com/technical_references/pdfs/electronic_engineering/Antenna_arrays.pdf) some intro notes to antenna arrays)

Types of arrays, depending on their current distribution:
- [Uniform array](uniform.m/): currents are the same in all elements of the array. Gives secondary lobes whenever the elements are more than half-wavelength apart.
- [Binomial array](binomial.m/): currents follow a binomial distribution. Gives no secondary lobes, only a wide primary lobe. 
- [Triangular array](triangular.m/): current ratios are triangular. Gives narrower main lobes than binomial, but with some secondary lobes (althought he secondary lobes are smaller than in the uniform array). 

[This file](Comparison.m/) compares all 3 current distributions mentioned above. 

It is possible to *move* the null points in the radiation pattern, see [this file](null_spacing.m/).

Finally [this file](null_spacing_Comparison.m/) makes a comprehensive comparison of any type of array mentioned above.
