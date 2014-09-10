This file contains miscellaneous notes about data processing.



Identifying missing output files
--------------------------------

```
from glob import glob
from os.path import join, exists

dirs = glob('2014.09.03-tax-parameter-sweep/mock-community/*/*/*/*/')

for e in dirs:
    fp = join(e,'rep_set_tax_assignments.txt')
    if not exists(fp):
        print e
```
