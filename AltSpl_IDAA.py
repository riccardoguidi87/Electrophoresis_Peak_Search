#### Python Script for the identification of peaks in IDAA assasy to be used as Alt-Splicing readout

from __future__ import division
import re
import os
import time
import math 
import string

## FUNCITONS ##
def Find_Peak(list1,wt_peak_size,tollerance):
#list1 should have the structure of the Dictionary value: |size,height|size,height|size,height|...
#wt_peak_size is inputed as a number for the search of the interested peak to be found in the Dictionary Value
#tollerance refers to the tollerance given in bp to the location of the peak of interest
#FUNCITON RETURNS the size,height of the peak seached within tollerance.
#V2.0 the function takes into considerastion that more then one peak can be found within the tollerance range. In that case, the
# program returns two potential peaks: one being the HIGHEST peak in the tollerange range, one being the closest to the
# exact SIZE of the expected peak
	
	list2 = []
	f = list1.strip("|").split("|")
	upper = wt_peak_size + tollerance
	lower = wt_peak_size - tollerance
	peak_search_1_size = 'not found'
	peak_search_1_height = 1	
	peak_search_2_size = 'not found'
	peak_search_2_height = 1
	distance = tollerance

	for i in f:
		# is the peak in range? if so, store height and size
		if  lower < float(i.split(",")[0]) < upper:
			list2.append(i)

	# with a list that contains all possible peaks do two searches:
	# for the HIGHEST PEAK in the range 
	for i in list2:
		if float(i.split(",")[1]) > peak_search_1_height:
			peak_search_1_size = float(i.split(",")[0])
			peak_search_1_height = float(i.split(",")[1])

	# for the closest to expected peak SIZE in the range
	for i in list2:
		new_distance = abs(float(i.split(",")[0]) - wt_peak_size)
		if new_distance < distance:
			distance = new_distance
			peak_search_2_size = float(i.split(",")[0])
			peak_search_2_height = float(i.split(",")[1])

	return str(peak_search_1_size) + "," + str(peak_search_1_height) + "," +  str(peak_search_2_size) + "," + str(peak_search_2_height)

def Find_Highest_Peak(list1):
#list1 should have the structure of the Dictionary value: |size,height|size,height|size,height|...
# FUNCTION RETURN the size of the highest peak in a PeakDict value: |size,height|size,height|size,height|...
# FUNCTION RETUNS 1 is no peaks found at all

	f = list1.strip("|").split("|")
	peak_size = 1
	peak_height = 1
	for i in f:
		if float(i.split(",")[1]) > peak_height:
			peak_size = float(i.split(",")[0])
			peak_height = float(i.split(",")[1])
	return peak_size

def Find_Peak_Pair(list1):
#list1 should have the structure of the Dictionary value: |size,height|size,height|size,height|...
#FUNCTION RETURNS two values, one for each of the two highest peaks in the run

	f = list1.strip("|").split("|")
	peak_size1 = 1
	peak_height1 = 1
	peak_size2 = 1
	peak_height2 = 1
	for i in f:
		if float(i.split(",")[1]) > peak_height1:
			peak_size1 = float(i.split(",")[0])
			peak_height1 = float(i.split(",")[1])
	
	for i in f:
		if float(i.split(",")[0]) != float(peak_size1):
			if float(i.split(",")[1]) > peak_height2 and float(i.split(",")[1]) < peak_height1 and abs(float(i.split(",")[0]) - float(peak_size1)) > 10:
				peak_size2 = float(i.split(",")[0])
				peak_height2 = float(i.split(",")[1])

	return str(peak_size1) + "," + str(peak_height1) + "," + str(peak_size2) + "," + str(peak_height2)


def KO_efficiency(list1,wt_peak_size):
#list1 should have the structure of the Dictionary value: |size,height|size,height|size,height|...
#wt_peak_size is the value obtained by Find_WT_Peak funciton
#return the relative abundance of the WT peak in the area surrounding the peak. This funciton requires 
#the identificaiton of the exact WT peak based on the WT run - see Find_WT_peak function
# note that the peak-range here is more restrictive then in the WT search, as we expect CRISPR cuts to be below 10bp

	f = list1.strip("|").split("|")
	upper = wt_peak_size + 10
	lower = wt_peak_size - 10
	sum_of_peaks = 1
	size = 1
	delta_in_hold = 100
	height = 1
	possible_wt_size = 1
	for i in f:
		# is the peak in range? if so, store height and store height in SUM of heights
		if  lower < float(i.split(",")[0]) < upper:
			sum_of_peaks += float(i.split(",")[1])
			delta = float(i.split(",")[0]) - float(wt_peak_size)
			delta = abs(delta)
			if delta < delta_in_hold: # here we test if the peak is closer to the WT size than a previously peak
				delta_in_hold = delta
				possible_wt_size = float(i.split(",")[0])
				height = float(i.split(",")[1])
	#if float(height) == 1:
		#print "although I found the WT peak in the WT sample, I didn't find the same peak in the KO sample, so I cannot calculate the relative KO efficienty"
	relative_abundance = float(height) / float(sum_of_peaks)
	return possible_wt_size, round(relative_abundance, 3)

##### END OF FUNCTIONS #####

# get info on how data are organized in the original Thermo table
with open('Table1.csv') as f: 
	line = f.readline().rstrip("\r\n").split(",")
	position_SampleName = line.index('"Sample File Name"')
	position_DyeColour = line .index('"Dye Color"')
	position_Size = line.index('"Size"')
	Position_Height = line.index('"Height"')

## create empty dictionary
Peaks_Dict = {}

## pre-populate Dict with all the Sample name files. This way you will create unique Key entries, AND add a value 1,1 to bypass an error
## Keep the Original Full Sample Name that include both info on Plate number and Well (necessary when data is a combo of multiple plate runs)

Table1 = open("Table1.csv")
next(Table1) # this comand skip the first line of the file
for line in Table1:
	SampleName = line.rstrip("\r\n").replace('"',"").split(",")[0]
	Peaks_Dict[SampleName] = "|1,1"
Table1.close()

## populate Dicrionary Peaks_Dict directly from the orignal Thermo table
Table1 = open ("Table1.csv")
next(Table1) # this comand skip the first line of the file
for line in Table1:
	linelist = line.rstrip("\r\n").replace('"',"").split(",")
	if linelist[position_DyeColour] == 'BLUE':
		if float((linelist[position_Size])) >= 100:
			if float((linelist[Position_Height])) >= 10:
				SampleName = linelist[0]
				Peaks_Dict[SampleName] = Peaks_Dict[SampleName] + "|" + str(linelist[position_Size]) + "," + str(linelist[Position_Height].rstrip('\n'))
Table1.close()

## at this stage, you should have Dictionary that contains one key per well per plate. Each Key has a value that correspond
## to the peak size and height (anything larger then 100bp and highr then 100) 
## Key = | size, height| size, height| size, height| size, height 

#### PART 2: Interrogate the Dictionary via iterations
## use a User_input table from which outsourse all the necessary info for the search

output = open("output.csv","w")
with open('User_Input.csv') as f:
    header = f.readline().rstrip("\r\n")
output.write(header)
output.write(",well#,PotentialPeak1_Search1_S,PotentialPeak1_Search1_H,PotentialPeak1_Search2_S,PotentialPeak1_Search2_H,PotentialPeak2_Search1_S,PotentialPeak2_Search1_H,PotentialPeak2_Search2_S,PotentialPeak2_Search2_H,HighestPeak_S,HighestPeak_S,HighestPeak_H,SecondHighestPeak_S,SecondHighestPeak_H,TotNumbPeaks\n") # estra-header for Question 1

## build a reference list that contains ordered wells in a 96WP format
## the list is [A1,A2,A3,A4... etc to H12]. The list is ordered.
## beware of the FORMAT OF THE WELLS NAMES: it is not 'A01' but 'A1'. This will change its usage!
## this requires the imprt string to work
n = 0
m = 0
plate96WP = []
for i in string.ascii_uppercase:
	row = str(i)
	for l in range(1,13):
		n += 1
		m += 1
		column = str(l)
		plate96WP.append(row + column)
		if n == 12:
			break
	if m == 96:
		break


User_input = open("User_Input.csv")
next(User_input)
for line in User_input:
	PlateID = line.split(",")[0]
	PPID = line.split(",")[1]
	expectedPeak1 = float(line.split(",")[2])
	expectedPeak2 = float(line.split(",")[3])
	Tollerance = float(line.split(",")[4])
	Gene = line.split(",")[5]
	FirstWell = line.split(",")[6]
	LastWell = line.split(",")[7].strip('\r\n')

	# which wells are interested 
	firstwell_location = plate96WP.index(FirstWell)
	lastwell_location = plate96WP.index(LastWell) + 1
	wells_range = plate96WP[firstwell_location:lastwell_location]
	
	# run loop
	for key, value in Peaks_Dict.items():
		if key.split("_")[0] == PlateID and key.split("_")[2] in wells_range:  ## CAREFUL HERE check SampleNames format!!
		#if key.split("_")[3] in wells_range:  ## USE THIS IF you have only one plate
			WellID = key.split("_")[2]
			Peak1 = Find_Peak(value,expectedPeak1,Tollerance)   ## Q1: Do I find anywhere (in the relevant plate and wells) the expected peaks?
			Peak2 = Find_Peak(value,expectedPeak2,Tollerance) 	## Q1: Do I find anywhere (in the relevant plate and wells) the expected peaks?
			HighestPeak = Find_Highest_Peak(value)
			Tot_Numb_Peaks = len(value.split("|")) - 2
			Peak_Pair = Find_Peak_Pair(value)
		
			# print "Looking in Plate ID: " + PlateID + " and looking in wellID: " + WellID

			output.write(str(PlateID) +','+ 
				str(PPID) +','+ 
				str(expectedPeak1) +','+ 
				str(expectedPeak2) +','+
				str(Tollerance) +','+
				str(Gene) +','+
				str(FirstWell) +','+
				str(LastWell) +','+
				str(WellID) +','+
				str(Peak1) +','+
				str(Peak2) + ',' +
				str(HighestPeak) + ',' +
				str(Peak_Pair) + ',' +
				str(Tot_Numb_Peaks) + '\n')

User_input.close()
output.close()
