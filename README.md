## Electrophoresis_Peak_Search
### A python script that search for desired peaks in an elecrophoresis table, to return peak height 

This program takes two tables as input:
1. **ThermoFisher Cloud Electrophoresis Tabulation** (TFET): you can obtain such table when you run your electrophoresis files (.fsa) using a [Cloud Thermo](https://www.thermofisher.com/order/catalog/product/A26811?SID=srch-srp-A26811#/A26811?SID=srch-srp-A26811) account, or any registered Thermo Electrophoresis softwares that can process .fsa file formats (that is the raw data from electrophoresis)
2. **User-Input-Table**: containing the necessary instructions for the identification of a specific peak, or set of peaks, and its quantification

Both tables are *.csv*

In the basic version of the program, the code searches for one peak, indicated in the User-Input-Table, inside the TFET. 

This program assumes you are running a 96WP-format of electrophoresis. You don't need to have a TFET containing all 96 sample, 
but each samples must be labelled with a row (A-H) and column system (1-12) that corresponds to a 96WP.

In the .fsa file, metadata contain the name of each sample, which will be passed to the TFET to label each individual sample. In this  version of the program, the nomenclature of each sample **must follow a precise order** (see below). The nomenclature I use for all my electrophoresis runs - and therefore the one this program can handle - looks like this:

```SampleName``` : "RUN200225B_A01_A1_015_extrainfo.fsa"

Where we have in order:

```SampleName``` :  *plateID _ wellcode _wellcode_samplenumber_extrainfo.fsa*

I am not sure why, but all the TFET I  contain a repetition of the ```wellcode``` (row and Col values combined). The program needs these two ```wellcode``` to be present in the ```SampleName``` in this peocire order.

If the ```SampleName``` is prepared differently, and cannot be easily changed, take a look at ```#run loop``` part of the Python script to fix the problem: 

```python
# run loop
	for key, value in Peaks_Dict.items():
		if key.split("_")[0] == PlateID and key.split("_")[2] in wells_range:  ## CAREFUL HERE check SampleNames format!!
		#if key.split("_")[3] in wells_range:  ## USE THIS IF you have only one plate
			WellID = key.split("_")[2]
```

... the location ```key.split("_")[2]``` of the split string refers to the location in ```SampleName``` where to find the ```wellcode``` in format ND (see above). If ```wellcode``` is present in another location, adjust accordingly. Also, if the split-value isn't "_", that can be changed here. Also ```key.split("__")[4]``` make sure this location in the SampleName hold the well ID in a format that is ‘A1’, not ‘A01’

WellID = key.split("_")[4]
make sure this IS the location of the well ID in SampleName a format that is ‘A1’, not ‘A01’

- The program does take into consideration you may have more then one plate in Table1.csv, and each plate may contain more then one PP. Make sure that in User_Input.csv you have correctly placed the plate ID (it does’t have to be a number, it can also be a string - most recent agreement in plate labeling from Dixon, the PlatesID are continuous letters A,B,C…. etc… ).

### Instructions

- A copy of the original TFET and the User-Input-Table must be present inside the same folder where the Python program runs
- Rename TFET as: “Table1.csv” - be exact in the name!
- The User-Input-Table must be named as "User_Input.csv", and should have this outlook:

### User-Input-Table Outlook

| PlateID    | PP#  | ExpectedPeak1 (nt) | ExpectedPeak2 (nt) | Tollerance (nt) | Gene  | FirstWell | LastWell |
| ---------- | ---- | ------------------ | ------------------ | --------------- | ----- | --------- | -------- |
| RUN200916A | A05  | 133                | 195                | 20              | Haus7 | A1        | H12      |
| RUN200916B | A07  | 128                | 164                | 20              | Ggta1 | A1        | H12      |
| RUN200916C | B01  | 253                | 417                | 20              | Tyk2  | A1        | H12      |

Notes on the User Input Table:

- Even if you have run only one electrophoresis plate (in the example above I run three), the name of plate must always be present in User Input Table, as well as in the ```SampleName``` of the TFET

- When modify the User_Input_Table using Excel, make sure you ERASE ANYTHING (even a blank cell) from underneath the last row (the program sometimes see the leftovers and tries to continue reading rows). Save table as *.csv*
- The program takes into consideration continues and ordered use of wells in the 96WP sent for electrophoresis - the program may still be working fine if some wells are to be missing from the run (as in, Table1.csv does not contains a raw for each of the 96 wells). 
- Beware of the way the facility lab that runs the IDAA plate labels the samples (usually first column in the exported Thermo Table): In the #run loop section, the program relies on a specific order of information present in the “Sample File Name” column of Table1.csv.



### Explanations of output table

well#
PotentialPeak1_Search1_S / *H = using first peak as a guidance, this is the peak S and H of the HIGHEST peak in tollerance
PotentialPeak1_Search2_S / *H = using first peak as guidance, this is the peak S and H with the SIZE closest to the expected size (yet this may not be the highest product)
PotentialPeak2_ect = as above, but for the second Alt-spl peak

HighestPeak_S, HighestPeak_H = the absolute highest peak in reaction (that passed Dictionary filtering) S and H

SecondHighestPeak_S,SecondHighestPeak_H = self explanatory

TotNumbPeaks = how many peaks the PCR has



















