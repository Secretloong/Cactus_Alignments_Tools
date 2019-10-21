#!/hwfssz5/ST_DIVERSITY/B10K/PUB/local/compilation_tool/Python-2.7.15/bin/python
import sys,fileinput

Argument = sys.argv[1:]
bed_file = Argument[0]
window = Argument[1]

#print bed_file

fd = fileinput.input(bed_file)

for line in fd:
	fa = line.rstrip("\n")
	end = int(fa.split()[2]) + 1
	for i in range(0,end,int(window)):
#		print str(end)+"\t"+str(i)
		if i+int(window) < int(fa.split()[2]):
			print fa.split()[0]+"\t"+str(i)+"\t"+str(i+int(window))+"\t"+fa.split()[0]+"_"+str(i)+"_"+str(i+int(window))
		elif i+int(window) == int(fa.split()[2]):
			print fa.split()[0]+"\t"+str(i)+"\t"+fa.split()[2]+"\t"+fa.split()[0]+"_"+str(i)+"_"+fa.split()[2]
		elif i+int(window) > int(fa.split()[2]) and i < int(fa.split()[2]):
			print fa.split()[0]+"\t"+str(i)+"\t"+fa.split()[2]+"\t"+fa.split()[0]+"_"+str(i)+"_"+fa.split()[2]
