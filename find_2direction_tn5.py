#!/usr/bin/env python

import sys
import regex
import pprint

sys.argv[1]
file_path = sys.argv[1]
file_handle = open(file_path, "r")
text = file_handle.read()

#define pattern match/count function
def count_matches(pattern):
	multiline_regex_object = regex.compile(pattern, regex.MULTILINE)
	count = 0
	for result in multiline_regex_object.finditer(text, overlapped=True):
		count += 1
	print(count)

ISpattern = r'\s[1-9]'
seqtwodirection = r"(^)([1-9][0-9]*?)( [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* )([1-9][0-9]*?)(\n)"
seqtwodirection2 = regex.compile(r"(^)([1-9][0-9]*?)( [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* )([1-9][0-9]*?)(\n)", regex.MULTILINE)
sixbpapart = regex.compile(r'^[1-9][0-9]*? [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [1-9][0-9]*?\n', regex.MULTILINE)
twelvebpapart = regex.compile(r'^[1-9][0-9]*? [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [0-9]*\n[0-9]* [1-9][0-9]*?\n', regex.MULTILINE)

if text.count("\n") < 50:
	print("Original plot file:	", text)
else:
	print("Plot file is too long to print:	", text.count("\n"), "lines.")
	
print("Total number of insertion sites: ",
count_matches(ISpattern))

print("Number of IS pairs 9bp apart, opposite directions (likely 2-direction sequencing of a single IS): ", 
count_matches(seqtwodirection))

print("Expected number based on insertion density:	")
sixbp = len(regex.findall(sixbpapart, text, overlapped=True))
twelvebp = len(regex.findall(twelvebpapart, text, overlapped=True))
expectednumber = (sixbp + twelvebp)/2
actualnumber = len(regex.findall(seqtwodirection2, text, overlapped=True))
print(expectednumber)
print(expectednumber/actualnumber)

def repl(m):
	first = str(int(m.group(2)) + int(m.group(4)))
	return m.group(1) + first + m.group(3) + str(0) + m.group(5)

if actualnumber/expectednumber > 5:
	print("High proportion of seq artefacts, collapsing to single site...")
	processed = regex.sub(seqtwodirection2, repl, text, overlapped=True)

	print("Total IS number after collapsing 2-direction sequencing artefacts:")
	newISnumber = regex.findall(ISpattern, processed)
	print(len(newISnumber))

	my_file = open(file_path + "_processed.txt", "w")
	print("Writing processed plot sites to file:	", file_path + "_processed.txt")
	my_file.write(processed)
	my_file.close()

else:
	print("Proportion less than 5:1, not going to collapse sites")
