'''
Created on 9 Jul 2013

@author: mcbryan
'''

import gzip
import mimetypes
import re

class FastQEntry(object):
    def __init__(self, header, sequence, qheader, quality):
        self.header = header
        self.sequence = sequence
        self.qheader = qheader
        self.quality = quality
        
        assert header.startswith("@")
        assert qheader.startswith("+")
    
    def minQual(self):
        return min([ord(base) for base in self.quality])
    
    def maxQual(self): 
        return max([ord(base) for base in self.quality])
    
    def headerFormat(self):
        
        # @HWUSI-EAS100R:6:73:941:1973#0/1
        # @<instrument>:<lane>:<tile>:<x>:<y>#<index>/<pair>
        
        if re.match("@.*:[1-8]:[0-9]*:[0-9]*:[0-9]*#.*/[12]",self.header):
            return "Illumina"
        
        # @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
        # @<instrument>:<runid>:<flowcell>:<lane>:<tile>:<x>:<y> <pair>:<filter>:<controlbits>:<index>
        if re.match("@.*:[0-9]*:.*:[1-8]:[0-9]*:[0-9]*:[0-9]* [12]:[YN]:[0-9]*:[GATC]*",self.header):
            return "Illumina1.8"
        
        return "Unknown"
        
    
    def __str__(self):
        return self.header + "\n" + self.sequence + "\n" + self.quality
        
class FastQFile(object):

    def __init__(self, filename):
        self.filename = filename
        
    def readlines(self):
        if self.filename.endswith("gz") or "gzip" in mimetypes.guess_type(self.filename):
            self.reader = gzip.open(self.filename)
        else:
            self.reader = open(self.filename)
        
        while True:
            
            four_lines = [self.reader.readline() for _ in range(4)]
            
            assert len(four_lines) == 4
            assert "\n" not in four_lines # we don't handle blank lines at the moment

            if "" in four_lines: break  # end of file
            
            # strip whitespace
            header,sequence,qheader,quality = [line.strip() for line in four_lines]
            
            yield FastQEntry(header, sequence, qheader, quality)
            
        self.reader.close()
    
    def __iter__(self):
        return self.readlines()        
        
if __name__ == "__main__":
    
    for i in FastQFile("/mnt/50tb/privatedata/David/MoreHistoneChIP/sequences/PD32.H4.lane1.fastq.gz"):
        print str(i)