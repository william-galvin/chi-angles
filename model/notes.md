# Notes

- Lowering batch size reduces overfitting, but does not 
improve validation scores

- Things I have tried varying: 
	- Number of CG blocks
	- ch_size_list
	- lmax_list
	- Learning rate
	- Batch size

- Loss Functions
	- Added noise to the invalid entries
	- Divisor term is number of valid entries

- Things to try next:
	- More data? (Twice as much? -- is there a new pipeline for this?)
