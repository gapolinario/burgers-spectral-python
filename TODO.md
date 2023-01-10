1. Burgers Initial, asserts were disabled. Why aren't they working? Should `tol` be bigger?
2. ensemble index `R` should not be accounted in the hashlib, this is fine for the deterministic case
3. Burgers White Noise, zero mode is not zero, this conflicts with an assertion in line 354 of the Burgers code