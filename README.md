## Electrophoresis_Peak_Search
A python script that search for approriate peaks in an elecrophoresis table

This program takes two tables as input: 
- ThermoFisher Cloud Electrophoresis Tabulation - you opbtain such talbe when you run your electrophoresis files (.fsa) with
 [Cloud Thermo](https://www.thermofisher.com/order/catalog/product/A26811?SID=srch-srp-A26811#/A26811?SID=srch-srp-A26811)

- Import in this folder a copy of the original Therm Table
- Rename the Original Table “Table1.csv” - be exact in the name!!
- Modify the User_Input.csv table to incorporate all the assays you have, and the specific location of each well

OBS: you very commonly fill up this table incompletely that results in errors!!

- When modify the User_Inout table, make sure you ERASE ANYTHING (even a blank cell) from underneath the last row (the program sometimes see the leftovers and tries to continue reading rows)
- The program takes into consideration continues and ordered use of wells in the 96WP sent for IDAA - the program may still be working fine if some wells are to be missing from the run (as in, Table1.csv does not contains a raw for each of the 96 wells). 
- Beware of the way the facility lab that runs the IDAA plate labels the samples (usually first column in the exported Thermo Table): In the #run loop section, the program relies on a specific order of information present in the “Sample File Name” column of Table1.csv. 
In particular, the program expect this:

RICC_B05_B5_045_2020-01-21.fsa
DATE_PLATE#_WELL#_ etc etc

Standardize labeling:
RUN_200225_B_A01_A1_015_date.fsa
RUN_plateID_plate.letter_Well_Well_number_YYYY-MM-DD.fsa

If the Sample File Name is named differently, take a look at #run loop part to fix it: 

key.split(“_”)[2]
make sure Sample Name is splittable using “_”, and that location [n] holds the number/letter of the plate

key.split(“_”)[4]
make sure this location in the SampleName hold the well ID in a format that is ‘A1’, not ‘A01’

WellID = key.split("_")[4]
make sure this IS the location of the well ID in SampleName a format that is ‘A1’, not ‘A01’

- The program does take into consideration you may have more then one plate in Table1.csv, and each plate may contain more then one PP. Make sure that in User_Input.csv you have correctly placed the plate ID (it does’t have to be a number, it can also be a string - most recent agreement in plate labeling from Dixon, the PlatesID are continuous letters A,B,C…. etc… ).

Explanations of output table
well#
PotentialPeak1_Search1_S / *H = using first peak as a guidance, this is the peak S and H of the HIGHEST peak in tollerance
PotentialPeak1_Search2_S / *H = using first peak as guidance, this is the peak S and H with the SIZE closest to the expected size (yet this may not be the highest product)
PotentialPeak2_ect = as above, but for the second Alt-spl peak

HighestPeak_S, HighestPeak_H = the absolute highest peak in reaction (that passed Dictionary filtering) S and H

SecondHighestPeak_S,SecondHighestPeak_H = self explanatory

TotNumbPeaks = how many peaks the PCR has















