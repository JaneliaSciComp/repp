1.2.2
- Made some changes in the enzyme digest - there were some inconsistencies in how a circular backbone was circularized
- Made enzymes optional - if no enzyme specified we lookup the backbone in the target
1.2.1
- The fragment IDs in the CSV output, are based on the input file instead of the target sequence ID
- Output filename is optional now - if not provided it is generated from the input
- Added support for CSV output (default output now). The CSV output consists of 2 files one for the cloning strategy and one for the reagents. The reagents can be populated from an input CSV file
- Add database supports adding files from one or multiple folders; errors will not
stop the program and they will be simply be reported at the end.