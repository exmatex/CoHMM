import subprocess
import re
import string

proc = subprocess.Popen(["ms_print", "massifOut"], stdout=subprocess.PIPE)
msOut = proc.stdout.read()
m = re.search('\d\s\(peak\)', msOut)
peakVal = m.group(0)
m = re.search('\d+', peakVal)
lineN = m.group(0)
regLine = r'\s+' + str(lineN) + r'\s+\d*,*\d*,*\d*,*\s+\d*,*\d*,*\d*,*\d*'
m = re.search(regLine, msOut)
goodLine = m.group(0)
splitLine = goodLine.split()
print splitLine[2]
